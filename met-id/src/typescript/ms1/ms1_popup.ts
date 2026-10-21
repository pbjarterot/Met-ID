import { invoke } from '@tauri-apps/api/core';
import { open } from '@tauri-apps/plugin-dialog';
import { createTableRowHTML } from './ms1_table';

let deleted_ms1_rows: number[] = [];
let dataset_transposed: boolean = false;
let glbl_result: string | string[] | null = "";
let inputType: string = "";

let ionMobilityPreviousState: {
    checked: boolean;
    disabled: boolean;
} | null = null;

let currentImzmlPath: string | null = null;
let imzmlUpdateTimer: number | undefined;
let imzmlRequestId = 0;

window.addEventListener("DOMContentLoaded", () => {
    document
        .getElementById("imzml-peak-count")
        ?.addEventListener("input", scheduleImzmlUpdate);

    document
        .getElementById("imzml-tolerance-ppm")
        ?.addEventListener("input", scheduleImzmlUpdate);
});

function getImzmlOptions(): { n: number; tolerancePpm: number } {
    const n = Number(
        (document.getElementById("imzml-peak-count") as HTMLInputElement).value
    );
    const tolerancePpm = Number(
        (document.getElementById("imzml-tolerance-ppm") as HTMLInputElement).value
    );

    if (!Number.isInteger(n) || n < 1) {
        throw new Error("Number of peaks must be a positive whole number.");
    }

    if (!Number.isFinite(tolerancePpm) || tolerancePpm < 0) {
        throw new Error("Tolerance must be a non-negative number.");
    }

    return { n, tolerancePpm };
}

async function read_input_imzml(
    path: string,
    n: number,
    tolerancePpm: number
): Promise<string[][]> {
    return invoke<string[][]>("read_ms1_imzml", {
        path,
        n,
        tolerancePpm,
        transpose: false
    });
}

function scheduleImzmlUpdate(): void {
    if (!currentImzmlPath) return;

    if (imzmlUpdateTimer !== undefined) {
        window.clearTimeout(imzmlUpdateTimer);
    }

    imzmlUpdateTimer = window.setTimeout(() => {
        void refreshImzmlSpectrum();
    }, 250);
}

async function refreshImzmlSpectrum(): Promise<void> {
    if (!currentImzmlPath) return;

    const requestId = ++imzmlRequestId;
    setImzmlPlotLoading(true);

    try {
        const { n, tolerancePpm } = getImzmlOptions();
        const rows = await read_input_imzml(
            currentImzmlPath,
            n,
            tolerancePpm
        );

        if (requestId !== imzmlRequestId) return;

        add_to_popup_table(rows);
        drawImzmlSpectrum(rows, tolerancePpm);
    } catch (error) {
        console.error("Could not update the imzML spectrum:", error);
    } finally {
        if (requestId === imzmlRequestId) {
            setImzmlPlotLoading(false);
        }
    }
}

function setIonMobilityForImzml(isImzml: boolean): void {
    const toggle = document.getElementById(
        "toggle-ion-mobility"
    ) as HTMLInputElement | null;

    if (!toggle) return;

    if (isImzml) {
        if (ionMobilityPreviousState === null) {
            ionMobilityPreviousState = {
                checked: toggle.checked,
                disabled: toggle.disabled
            };
        }

        toggle.checked = false;
        toggle.disabled = true;
    } else {
        if (ionMobilityPreviousState !== null) {
            toggle.checked = ionMobilityPreviousState.checked;
            toggle.disabled = ionMobilityPreviousState.disabled;
            ionMobilityPreviousState = null;
        }
    }

    // Let the existing switch listener update column visibility.
    toggle.dispatchEvent(new Event("change", { bubbles: true }));
}

async function data(result: string, result_type: string) {
    if (!result || result.length === 0) return;

    if (deleted_ms1_rows.length > 0) {
        deleted_ms1_rows.length = 0;
    }

    const isImzml =
        result_type === "filename" &&
        result.toLowerCase().endsWith(".imzml");

    setIonMobilityForImzml(isImzml);

    if (isImzml) {
        currentImzmlPath = result;
    } else {
        currentImzmlPath = null;
    }

    // Ignore pasted content that does not look like delimited data.
    if (result_type !== "filename" && !containsDelimiter(result)) {
        return;
    }

    const csvControls = document.getElementById("ms1-popup-csv-controls");
    const imzmlControls = document.getElementById("ms1-popup-imzml-controls");

    if (csvControls) csvControls.hidden = isImzml;
    if (imzmlControls) imzmlControls.hidden = !isImzml;

    showPopup();
    setMs1Loading(true);

    try {
        let imzmlTolerancePpm = 10;
        const delimiter = ";";
        let fileContents: string[][] = [[]];

        if (result_type === "filename") {
            inputType = "file";

            if (isImzml) {
                const { n, tolerancePpm } = getImzmlOptions();
                imzmlTolerancePpm = tolerancePpm;

                fileContents = await read_input_imzml(
                    result,
                    n,
                    tolerancePpm,
                );
            } else {
                fileContents = await read_input_csv(
                    result,
                    "ccs",
                    20,
                    delimiter,
                    dataset_transposed,
                    deleted_ms1_rows
                );
            }
        } else {
            inputType = "ctrlv";
            fileContents = await read_ctrl_v(
                result,
                delimiter,
                dataset_transposed,
                deleted_ms1_rows
            );
        }

        add_filename_to_top_of_popup(result);
        add_to_popup_table(fileContents);

        if (isImzml) {
            drawImzmlSpectrum(fileContents, imzmlTolerancePpm);
        }


        // These are CSV column controls; leave them out for imzML.
        if (!isImzml) {
            addPopupDropdownElements();
        }

        glbl_result = result;
    } catch (error) {
        console.error("Failed to load MS1 data:", error);
        alert(`Failed to load data: ${String(error)}`);
    } finally {
        setMs1Loading(false);
    }
}

function setImzmlPlotLoading(loading: boolean): void {
    const spinner = document.getElementById("imzml-plot-loading");
    const chart = document.getElementById("ms1-imzml-chart-container");

    if (spinner) spinner.hidden = !loading;
    chart?.setAttribute("aria-busy", String(loading));
}

function setMs1Loading(loading: boolean): void {
    const spinner = document.getElementById("ms1-popup-loading");
    if (spinner) spinner.hidden = !loading;

    const continueButton = document.getElementById(
        "ms1-popup-continue-button"
    ) as HTMLButtonElement | null;

    if (continueButton) continueButton.disabled = loading;
}

export function get_ctrl_v_data(result:string) {
    data(result, "content");
}

export async function get_csv() {
    let result: string = await open({directory: false, multiple: false, 
        filters: [{name: 'Supported Files', extensions:['csv', 'tsv', 'imzml', 'txt']}]}) as string;
    data(result, "filename");
    
}

function containsDelimiter(str: string): boolean {
    return str.indexOf(',') !== -1 || str.indexOf(';') !== -1 || str.indexOf('\t') !== -1;
}

async function read_input_csv(result: string, typeofcsv:string, n: number, delimiter: string, transpose: boolean, deleted_rows: number[]) {
    const parsed_csv: string[][] = await invoke("read_ms1_csv", {path: result, typeofcsv:typeofcsv, n:n, delimiter: delimiter, transpose: transpose, deletedrows: deleted_rows});
    console.log("parsed_csv", n, parsed_csv)
    return parsed_csv
}
/*
async function read_input_imzml(result: string, transpose: boolean) {
    const fileContents: string[][] = await invoke("read_ms1_imzml", {path: result, n: 1000, tolerancePpm: 10, transpose: transpose});
    drawImzmlSpectrum(fileContents)
    return fileContents
}
*/
async function read_ctrl_v(result: string, delimiter: string, transpose: boolean, deleted_rows: number[]) {
    const parsed_csv: string[][] = await invoke("read_ms1_ctrlv", {content: result, delimiter: delimiter, transpose: transpose, deletedrows: deleted_rows});
    return parsed_csv
}

import * as d3 from "d3";

interface Peak {
    mz: number;
    intensity: number;
}

let spectrumResizeObserver: ResizeObserver | null = null;


function drawImzmlSpectrum(rows: string[][], tolerancePpm: number): void {
    const container = document.getElementById(
        "ms1-imzml-chart-container"
    );
    const svgElement = document.getElementById(
        "ms1-imzml-spectrum"
    ) as SVGSVGElement | null;

    if (!container || !svgElement) {
        console.error("Could not find the imzML spectrum container or SVG");
        return;
    }

    const peaks: Peak[] = rows
        .map(([mz, intensity]) => ({
            mz: Number(mz),
            intensity: Number(intensity)
        }))
        .filter(
            peak =>
                Number.isFinite(peak.mz) &&
                Number.isFinite(peak.intensity)
        );

    if (peaks.length === 0) {
        console.warn("No valid spectrum peaks to draw");
        return;
    }

    spectrumResizeObserver?.disconnect();

    const svg = d3.select(svgElement);

    function render(): void {
        const { width, height } = container!.getBoundingClientRect();

        // The popup may not have finished laying out yet.
        if (width <= 0 || height <= 0) return;

        const margin = { top: 20, right: 20, bottom: 55, left: 70 };
        const plotWidth = width - margin.left - margin.right;
        const plotHeight = height - margin.top - margin.bottom;

        if (plotWidth <= 0 || plotHeight <= 0) return;

        const minMz = d3.min(peaks, peak => peak.mz) ?? 0;
        const maxMz = d3.max(peaks, peak => peak.mz) ?? 1;
        const fullMaxIntensity =
            d3.max(peaks, peak => peak.intensity) ?? 1;

        const x = d3.scaleLinear()
            .domain(
                minMz === maxMz
                    ? [minMz - 1, maxMz + 1]
                    : [minMz, maxMz]
            )
            .nice()
            .range([0, plotWidth]);

        // Intensity is recalculated against the tallest peak in the
        // currently visible m/z range.
        const y = d3.scaleLinear()
            .domain([0, 100])
            .range([plotHeight, 0]);

        svg
            .attr("width", width)
            .attr("height", height)
            .attr("viewBox", `0 0 ${width} ${height}`)
            .selectAll("*")
            .remove();

        const defs = svg.append("defs");

        const clipId = "ms1-imzml-spectrum-clip";

        defs.append("clipPath")
            .attr("id", clipId)
            .append("rect")
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", plotWidth)
            .attr("height", plotHeight);

        const plot = svg.append("g")
            .attr(
                "transform",
                `translate(${margin.left},${margin.top})`
            );

        const xAxis = plot.append("g")
            .attr("transform", `translate(0,${plotHeight})`);

        const yAxis = plot.append("g");

        const peaksGroup = plot.append("g")
            .attr("clip-path", `url(#${clipId})`);

        const pickedSpans = plot.append("g")
            .attr("clip-path", `url(#${clipId})`);

        peaksGroup.selectAll<SVGLineElement, Peak>("line")
            .data(peaks)
            .join("line")
            .attr("class", "peak")
            .attr("x1", peak => x(peak.mz))
            .attr("x2", peak => x(peak.mz))
            .attr("stroke", "#7db7ff")
            .attr("stroke-width", 2);

        plot.append("text")
            .attr("x", plotWidth / 2)
            .attr("y", plotHeight + 45)
            .attr("text-anchor", "middle")
            .text("m/z");

        plot.append("text")
            .attr("transform", "rotate(-90)")
            .attr("x", -plotHeight / 2)
            .attr("y", -50)
            .attr("text-anchor", "middle")
            .text("Relative intensity");

        function update(transform: d3.ZoomTransform): void {
            const zoomedX = transform.rescaleX(x);
            const [domainMin, domainMax] = zoomedX.domain();

            const visiblePeaks = peaks.filter(
                peak =>
                    peak.mz >= domainMin &&
                    peak.mz <= domainMax
            );

            const visibleMaxIntensity =
                d3.max(visiblePeaks, peak => peak.intensity) ??
                fullMaxIntensity;

            const intensityScale =
                visibleMaxIntensity > 0
                    ? visibleMaxIntensity
                    : 1;

            pickedSpans.selectAll<SVGRectElement, Peak>("rect")
                .data(visiblePeaks, peak => String(peak.mz))
                .join("rect")
                .attr("x", peak => {
                    const halfWidth = Math.abs(peak.mz) * tolerancePpm / 1_000_000;
                    return zoomedX(peak.mz - halfWidth);
                })
                .attr("y", plotHeight - 8)
                .attr("width", peak => {
                    const halfWidth = Math.abs(peak.mz) * tolerancePpm / 1_000_000;
                    const spanWidth =
                        zoomedX(peak.mz + halfWidth) -
                        zoomedX(peak.mz - halfWidth);

                    // Keep very narrow ppm spans visible on screen.
                    return Math.max(2, spanWidth);
                })
                .attr("height", 8)
                .attr("fill", "red")
                .attr("opacity", 0.8);

            xAxis.call(d3.axisBottom(zoomedX).ticks(Math.max(2, Math.floor(plotWidth / 50))));

            yAxis.call(
                d3.axisLeft(y)
                    .tickFormat(value => `${value}%`)
            );

            peaksGroup.selectAll<SVGLineElement, Peak>("line.peak")
                .attr("x1", peak => zoomedX(peak.mz))
                .attr("x2", peak => zoomedX(peak.mz))
                .attr(
                    "y1",
                    peak =>
                        y((peak.intensity / intensityScale) * 100)
                )
                .attr("y2", y(0));
        }

        // This transparent rectangle captures the wheel and drag events
        // over the plotting area.
        const zoomSurface = plot.append("rect")
            .attr("width", plotWidth)
            .attr("height", plotHeight)
            .attr("fill", "transparent")
            .style("pointer-events", "all");

        const zoom = d3.zoom<SVGRectElement, unknown>()
            .scaleExtent([1, 10000])
            .extent([
                [0, 0],
                [plotWidth, plotHeight]
            ])
            .filter(event =>
                event.type === "wheel" ||
                (event.type === "mousedown" && event.button === 0)
            )
            .on("zoom", event => {
                update(event.transform);
            });

        zoomSurface
            .call(zoom)
            // Remove D3's built-in double-click zoom.
            .on("dblclick.zoom", null)
            // Reset to the full spectrum on double-click.
            .on("dblclick.reset", event => {
                event.preventDefault();
                zoomSurface.call(zoom.transform, d3.zoomIdentity);
            });

        update(d3.zoomIdentity);
    }

    render();

    spectrumResizeObserver = new ResizeObserver(() => {
        render();
    });

    spectrumResizeObserver.observe(container);
}

export function showPopup() {
    const overlay = document.getElementById('overlay');
    const popup = document.getElementById('popup');

    overlay!.style.display = 'block';
    popup!.style.display = 'flex';
    
    const closePopupButton = document.getElementById("ms1-popup-cancel-button");
    closePopupButton?.addEventListener("click", () => hidePopup())

    const overlayDiv = document.getElementById("overlay");
    overlayDiv?.addEventListener("click", () => hidePopup());

    eventListenerModule.addPopupTransposeListener();
    eventListenerModule.addPopupChangeDelimiterListener();
    eventListenerModule.addPopupContinueListener();
    
    const delimiterInput = document.getElementById(
        "ms1-popup-delimiter-input"
    ) as HTMLInputElement | null;

    if (delimiterInput) {
        delimiterInput.value = ";";
    }
}

function hidePopup() {
    const overlay = document.getElementById('overlay');
    const popup = document.getElementById('popup');

    overlay!.style.display = 'none';
    popup!.style.display = 'none';

    eventListenerModule.removePopupTransposeListener();
    eventListenerModule.removePopupChangeDelimiterListener();
    eventListenerModule.removePopupContinueListener()
}

function add_popup_table_row(items: string[], i: number) {
    let template = `
    <tr id="ms1-popup-row${i}">
        <td><button class="ms1-popup-trash-button" id="ms1-popup-trash-button${i}"><ion-icon name="trash-outline"></ion-icon></button></td>
        <td>${items[0]}</td>
        <td>${items[1]}</td>
        <td>${items[2]}</td>
        <td>${items[3]}</td>
        <td>${items[4]}</td>
        <td>${items[5]}</td>
        <td>${items[6]}</td>
    </tr>
    `
    return template

}

type Ms1ColumnRole = "" | "m/z" | "CCS";

function addPopupDropdownElements(): void {
    const table = document.getElementById(
        "ms1-popup-datatable"
    ) as HTMLTableElement | null;

    const headerRow = table?.tHead?.rows[0] ?? table?.rows[0];
    if (!headerRow) return;

    const ionMobilitySwitch = document.getElementById(
        "toggle-ion-mobility"
    ) as HTMLInputElement | null;

    // Add this listener only once, even if the popup is opened repeatedly.
    if (
        ionMobilitySwitch &&
        ionMobilitySwitch.dataset.ms1PopupListener !== "true"
    ) {
        ionMobilitySwitch.addEventListener("change", updateMs1CcsOptions);
        ionMobilitySwitch.dataset.ms1PopupListener = "true";
    }

    const roles: Ms1ColumnRole[] = ["", "m/z", "CCS"];
    const usedRoles = new Set<Ms1ColumnRole>();

    // Cell 0 is the delete button; the remaining cells map to CSV columns.
    for (let cellIndex = 1; cellIndex < headerRow.cells.length; cellIndex++) {
        const cell = headerRow.cells[cellIndex];
        const columnIndex = cellIndex - 1;
        const previous = cell.querySelector<HTMLSelectElement>(
            ".ms1-popup-column-role"
        );

        let selectedRole = (previous?.value ?? "") as Ms1ColumnRole;

        if (!selectedRole) {
            if (columnIndex === 0) selectedRole = "m/z";
            if (columnIndex === 4 && ionMobilitySwitch?.checked) {
                selectedRole = "CCS";
            }
        }

        // Each role can be assigned to only one column.
        if (selectedRole && usedRoles.has(selectedRole)) {
            selectedRole = "";
        }
        if (selectedRole) usedRoles.add(selectedRole);

        const select = document.createElement("select");
        select.className = "ms1-popup-column-role";
        select.dataset.columnIndex = String(columnIndex);

        for (const role of roles) {
            const option = document.createElement("option");
            option.value = role;
            option.textContent = role || "—";
            select.appendChild(option);
        }

        select.value = selectedRole;

        select.addEventListener("change", () => {
            const role = select.value as Ms1ColumnRole;
            if (!role) return;

            document
                .querySelectorAll<HTMLSelectElement>(
                    ".ms1-popup-column-role"
                )
                .forEach(other => {
                    if (other !== select && other.value === role) {
                        other.value = "";
                    }
                });
        });

        cell.replaceChildren(select);
    }

    updateMs1CcsOptions();
}

function updateMs1CcsOptions(): void {
    const enabled = (
        document.getElementById("toggle-ion-mobility") as HTMLInputElement | null
    )?.checked ?? false;

    document
        .querySelectorAll<HTMLSelectElement>(".ms1-popup-column-role")
        .forEach(select => {
            const ccsOption = Array.from(select.options)
                .find(option => option.value === "CCS");

            if (ccsOption) ccsOption.disabled = !enabled;
            if (!enabled && select.value === "CCS") select.value = "";
        });
}

function getMs1PopupColumnMapping(): { mz: number; ccs?: number } {
    const mapping: { mz?: number; ccs?: number } = {};

    document
        .querySelectorAll<HTMLSelectElement>(".ms1-popup-column-role")
        .forEach(select => {
            const columnIndex = Number(select.dataset.columnIndex);

            if (select.value === "m/z") mapping.mz = columnIndex;
            if (select.value === "CCS") mapping.ccs = columnIndex;
        });

    if (mapping.mz === undefined) {
        throw new Error("Select one CSV column as m/z before continuing");
    }

    const useIonMobility = (
        document.getElementById("toggle-ion-mobility") as HTMLInputElement | null
    )?.checked ?? false;

    if (useIonMobility && mapping.ccs === undefined) {
        throw new Error("Select one CSV column as CCS before continuing");
    }

    return mapping as { mz: number; ccs?: number };
}

const popupTransposeButton = document.getElementById("ms1-popup-transpose-button");
async function popupTransposeListener(_event: Event) {
    dataset_transposed = !dataset_transposed;
    let delimiter = (document.getElementById("ms1-popup-delimiter-input") as HTMLInputElement).value as string;
    let file_contents: string[][] = [[]];
    if (inputType === "file") {
        file_contents = await read_input_csv(glbl_result as string, "ms1", 7, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    } else if (inputType === "ctrlv"){
        file_contents = await read_ctrl_v(glbl_result as string, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    }
    add_to_popup_table(file_contents)
}

const popupChangeDelimiterButton = document.getElementById("ms1-popup-change-delimiter-button");
async function popupChangeDelimiterListener(_event: Event) {
    let delimiter = (document.getElementById("ms1-popup-delimiter-input") as HTMLInputElement).value as string;
    let file_contents: string[][] = [[]];
    if (inputType === "file") {
        file_contents = await read_input_csv(glbl_result as string, "ms1", 7, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    } else if (inputType === "ctrlv"){
        file_contents = await read_ctrl_v(glbl_result as string, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    }
    
    add_to_popup_table(file_contents);
}

const popupContinueButton = document.getElementById("ms1-popup-continue-button");
async function popupContinueListener(_event: Event): Promise<void> {
    const delimiter = (
        document.getElementById("ms1-popup-delimiter-input") as HTMLInputElement
    ).value;

    try {
        await read_csv(glbl_result, deleted_ms1_rows, delimiter);
        hidePopup();
    } catch (error) {
        console.error("Failed to move popup data to the MS1 table:", error);
        alert(String(error));
    }
}

const eventListenerModule = {
    //Continue
    addPopupContinueListener() {
        popupContinueButton!.addEventListener("click", popupContinueListener);
    },
    removePopupContinueListener() {
        popupContinueButton!.addEventListener("click", popupContinueListener);
    },


    //Transpose
    addPopupTransposeListener() {
        popupTransposeButton!.addEventListener("click", popupTransposeListener);
    },
    removePopupTransposeListener() {
        popupTransposeButton!.removeEventListener("click", popupTransposeListener);
    },

    //Changing Delimiter
    addPopupChangeDelimiterListener() {
        popupChangeDelimiterButton!.addEventListener("click", popupChangeDelimiterListener)
    },

    removePopupChangeDelimiterListener() {
        popupChangeDelimiterButton!.removeEventListener("click", popupChangeDelimiterListener)
    },
};

function add_to_popup_table(file_contents) {
    const tbody = document.getElementById("ms1-popup-tbody")
        tbody!.innerHTML = ""
        let upper_limit = 50
        if (file_contents.length < 50) {
            upper_limit = file_contents.length
        }
        for (let i = 0; i < upper_limit; i++) {
            tbody!.insertAdjacentHTML("beforeend", add_popup_table_row(file_contents[i], i));
            let trashbutton = document.getElementById(`ms1-popup-trash-button${i}`)
            trashbutton!.addEventListener("click", delete_popup_table_row);
        }

}

function add_filename_to_top_of_popup(result){
    let resultString: string | null = Array.isArray(result) ? result.join(", ") : result;
    let filenameWindow = document.getElementById("ms1-popup-filename-window-name");
    filenameWindow!.textContent = resultString;
}

function delete_popup_table_row(event: Event) {
    const clickedButtonId = (event.target as HTMLElement).parentElement!.id;
    let row_idx = clickedButtonId.replace("ms1-popup-trash-button", "");

    const row = document.getElementById("ms1-popup-row" + row_idx as string);
    if (row) {
        row.remove();
    }
    deleted_ms1_rows.push(parseInt(row_idx));
}

function removeItemsByIndices(arr: number[], indicesToRemove: number[]): number[] {
    return arr.filter((_, index) => !indicesToRemove.includes(index));
}

async function read_csv(
    result: string | string[] | null,
    deleted_rows: number[],
    delimiter: string
): Promise<void> {
    const tableBody = document.getElementById("ms1-table-body");
    if (!tableBody) return;

    tableBody.innerHTML = "";

    const isImzml =
        inputType === "file" &&
        typeof result === "string" &&
        result.toLowerCase().endsWith(".imzml");

    if (isImzml) {
        const { n, tolerancePpm } = getImzmlOptions();
        const rows = await invoke<string[][]>("read_ms1_imzml", {
            n,
            tolerancePpm,
            path: result,
            transpose: false
        });

        const deleted = new Set(deleted_rows);
        let allRowsHTML = "";

        rows.forEach((row, index) => {
            if (deleted.has(index)) return;

            const mz = Number(row[0]);
            if (!Number.isFinite(mz)) return;

            // imzML provides m/z and intensity here, not CCS.
            allRowsHTML += createTableRowHTML(mz, Number.NaN, false);
        });

        tableBody.insertAdjacentHTML("beforeend", allRowsHTML);
        return;
    }

    // Only CSV and pasted data need the popup column mapping.
    const mapping = getMs1PopupColumnMapping();
    const useIonMobility =
        (document.getElementById("toggle-ion-mobility") as HTMLInputElement | null)
            ?.checked ?? false;

    let parsedCsv: number[][] = [[]];

    if (inputType === "file" && typeof result === "string") {
        parsedCsv = await invoke("read_ms1_csv", {
            path: result,
            typeofcsv: "ccs",
            n: 100000,
            delimiter,
            transpose: !dataset_transposed,
            deletedrows: []
        });
    } else if (inputType === "ctrlv") {
        parsedCsv = await invoke("read_ms1_ctrlv", {
            content: result,
            delimiter,
            transpose: !dataset_transposed,
            deletedrows: []
        });
    }

    const mzValues = parsedCsv[mapping.mz] ?? [];
    const ccsValues =
        useIonMobility && mapping.ccs !== undefined
            ? parsedCsv[mapping.ccs] ?? []
            : [];

    const masses = removeItemsByIndices(mzValues, deleted_rows);
    const ccss = useIonMobility
        ? removeItemsByIndices(ccsValues, deleted_rows)
        : [];

    const allRowsHTML = masses
        .map((mz, index) => {
            const ccs = useIonMobility
                ? ccss[index] ?? Number.NaN
                : Number.NaN;

            return createTableRowHTML(mz, ccs, useIonMobility);
        })
        .join("");

    tableBody.insertAdjacentHTML("beforeend", allRowsHTML);
}
