import { invoke } from '@tauri-apps/api/core';
import { fill_dropdown } from "../dropdown";
import { open} from '@tauri-apps/plugin-dialog';

//const dropdown_adducts = {"[M+H]+": 1, "[M+Na]+": 1, "[M+K]+": 1, "[M]+": 1, "[M-H]-": -1};
const dropdown_adducts_string = ["[M+H]+", "[M+Na]+", "[M+K]+", "[M]+", "[M-H]-"]
//let imageData = {"names": [[]] as string[][], "smiles": [[]] as string[][], "hmdb": [[]] as string[][]}

let ccsPredictionRunning = false;
let ccsPredictionPaused = false;
let ccsPredictionCancelled = false;

let ccsResumeResolver: (() => void) | null = null;

let ccsStartedAt = 0;
let ccsPausedAt: number | null = null;
let ccsTotalPausedTime = 0;
let ccsElapsedTimer: number | undefined;



let dataset_transposed: boolean = false;
//let inputType: string = "";
let deleted_ms1_rows: number[] = [];
let glbl_result: string | string[] | null = "";


window.addEventListener("DOMContentLoaded", async () => {
    //renderCCSHTML();
    setupStructureImageOverlay();

    fill_dropdown(dropdown_adducts_string, "dropdown-adducts");


    document.getElementById("predict-single-smiles-button")!
        .addEventListener("click", async () => ccs_prediction());

    document.getElementById("open-ccs-file-button")!
        .addEventListener("click", async () => open_ccs_file());

    document
        .getElementById("pause-ccs-prediction-button")
        ?.addEventListener("click", toggleCcsPredictionPause);

    document
        .getElementById("cancel-ccs-prediction-button")
        ?.addEventListener("click", cancelCcsPrediction);

    document.getElementById("save-ccs-prediction-button")!
        .addEventListener("click", saveCcsPredictions);

    document
        .getElementById("ccs-table-search")
        ?.addEventListener("input", filterCcsTable);

    const ccsSwitch = document.getElementById(
        "toggle-ccs-column"
    ) as HTMLInputElement | null;

    document
    .getElementById("clear-ccs-table-button")
    ?.addEventListener("click", clearCcsTable);

    ccsSwitch?.addEventListener("change", () => {
        toggleCcsColumn(ccsSwitch.checked);
    });

    if (ccsSwitch) {
        console.log("hi")
        toggleCcsColumn(ccsSwitch.checked);
    }
});

function clearCcsTable(): void {
    const tableBody = document.getElementById("ccs-table-body");

    if (!tableBody) {
        return;
    }

    tableBody.replaceChildren();

    // Also reset the search field.
    const searchInput = document.getElementById(
        "ccs-table-search"
    ) as HTMLInputElement | null;

    if (searchInput) {
        searchInput.value = "";
    }
}

function rowMatchesCcsSearch(
    row: HTMLTableRowElement,
    query: string
): boolean {
    const terms = query
        .trim()
        .toLowerCase()
        .split(/\s+/)
        .filter(Boolean);

    if (terms.length === 0) {
        return true;
    }

    // Only search the five text columns, not the structure preview.
    const rowText =
        row.dataset.searchText ??
        Array.from(row.cells)
            .slice(0, 5)
            .map(cell => cell.textContent ?? "")
            .join(" ")
            .toLowerCase();

    return terms.every(term => rowText.includes(term));
}

function filterCcsTable(): void {
    const input = document.getElementById(
        "ccs-table-search"
    ) as HTMLInputElement | null;

    const tableBody = document.getElementById("ccs-table-body");

    if (!input || !tableBody) {
        return;
    }

    tableBody
        .querySelectorAll<HTMLTableRowElement>("tr")
        .forEach(row => {
            row.style.display = rowMatchesCcsSearch(row, input.value)
                ? ""
                : "none";
        });
}

function formatElapsedTime(milliseconds: number): string {
    const totalSeconds = Math.floor(milliseconds / 1000);

    const hours = Math.floor(totalSeconds / 3600);
    const minutes = Math.floor((totalSeconds % 3600) / 60);
    const seconds = totalSeconds % 60;

    const paddedMinutes = String(minutes).padStart(2, "0");
    const paddedSeconds = String(seconds).padStart(2, "0");

    if (hours > 0) {
        return (
            `${String(hours).padStart(2, "0")}:` +
            `${paddedMinutes}:${paddedSeconds}`
        );
    }

    return `${paddedMinutes}:${paddedSeconds}`;
}

function updatePauseButton(paused: boolean): void {
    const button = document.getElementById(
        "pause-ccs-prediction-button"
    ) as HTMLButtonElement | null;

    if (!button) {
        return;
    }

    let icon = button.querySelector("ion-icon");

    if (!icon) {
        icon = document.createElement("ion-icon");
        button.replaceChildren(icon);
    }

    icon.setAttribute(
        "name",
        paused
            ? "play-circle-outline"
            : "pause-circle-outline"
    );

    const action = paused ? "Resume" : "Pause";

    button.setAttribute(
        "aria-label",
        `${action} CCS prediction`
    );

    button.title = `${action} CCS prediction`;
}

function updateCcsElapsedTime(): void {
    const elapsedElement =
        document.getElementById("ccs-elapsed-time");

    if (!elapsedElement || ccsStartedAt === 0) {
        return;
    }

    const now = performance.now();

    const currentPauseDuration =
        ccsPausedAt === null
            ? 0
            : now - ccsPausedAt;

    const elapsed =
        now -
        ccsStartedAt -
        ccsTotalPausedTime -
        currentPauseDuration;

    elapsedElement.textContent =
        `Elapsed: ${formatElapsedTime(Math.max(0, elapsed))}`;
}

function startCcsElapsedTimer(): void {
    if (ccsElapsedTimer !== undefined) {
        window.clearInterval(ccsElapsedTimer);
    }

    ccsStartedAt = performance.now();
    ccsPausedAt = null;
    ccsTotalPausedTime = 0;

    updateCcsElapsedTime();

    ccsElapsedTimer = window.setInterval(
        updateCcsElapsedTime,
        500
    );
}

function stopCcsElapsedTimer(): void {
    updateCcsElapsedTime();

    if (ccsElapsedTimer !== undefined) {
        window.clearInterval(ccsElapsedTimer);
        ccsElapsedTimer = undefined;
    }
}

function releaseCcsPause(): void {
    if (ccsPausedAt !== null) {
        ccsTotalPausedTime +=
            performance.now() - ccsPausedAt;
    }

    ccsPredictionPaused = false;
    ccsPausedAt = null;

    const resolver = ccsResumeResolver;
    ccsResumeResolver = null;
    resolver?.();
}

function toggleCcsPredictionPause(): void {
    if (!ccsPredictionRunning) {
        return;
    }

    const status = document.getElementById(
        "ccs-prediction-status"
    );

    if (!ccsPredictionPaused) {
        ccsPredictionPaused = true;
        ccsPausedAt = performance.now();

        updatePauseButton(true);

        if (status) {
            status.textContent = "Paused";
        }
    } else {
        releaseCcsPause();
        updatePauseButton(false);

        if (status) {
            status.textContent = "Running";
        }
    }
}

function cancelCcsPrediction(): void {
    if (!ccsPredictionRunning) {
        return;
    }

    ccsPredictionCancelled = true;

    // Release the loop if it is currently paused.
    if (ccsPredictionPaused) {
        releaseCcsPause();
    }

    const status = document.getElementById(
        "ccs-prediction-status"
    );

    if (status) {
        status.textContent = "Cancelling…";
    }
}

function waitWhileCcsPaused(): Promise<void> {
    if (!ccsPredictionPaused) {
        return Promise.resolve();
    }

    return new Promise<void>(resolve => {
        ccsResumeResolver = resolve;
    });
}


function setupStructureImageOverlay(): void {
    const table = document.getElementById("ccs-prediction-table");

    if (!table) {
        console.error("Could not find #ccs-prediction-table");
        return;
    }

    const overlay = document.createElement("div");
    overlay.id = "ccs-image-overlay";
    overlay.className = "ccs-image-overlay";
    overlay.hidden = true;

    overlay.innerHTML = `
        <div
            class="ccs-image-overlay-content"
            role="dialog"
            aria-modal="true"
            aria-label="Molecular structure"
        >
            <button
                type="button"
                class="ccs-image-overlay-close"
                aria-label="Close structure image"
            >
                &times;
            </button>

            <img
                id="ccs-image-overlay-image"
                alt="Enlarged molecular structure"
            >
        </div>
    `;

    document.body.appendChild(overlay);

    const enlargedImage = overlay.querySelector(
        "#ccs-image-overlay-image"
    ) as HTMLImageElement;

    const closeButton = overlay.querySelector(
        ".ccs-image-overlay-close"
    ) as HTMLButtonElement;

    function closeOverlay(): void {
        overlay.hidden = true;
        enlargedImage.removeAttribute("src");
        document.body.classList.remove("overlay-open");
    }

    table.addEventListener("click", event => {
        const target = event.target as HTMLElement;

        const clickedImage = target.closest<HTMLImageElement>(
            ".preview-image-cell img"
        );

        if (!clickedImage) {
            return;
        }

        const smiles = clickedImage.dataset.sourceSmiles;

        if (!smiles) {
            console.error("No source SMILES found for image");
            return;
        }

        overlay.hidden = false;
        document.body.classList.add("overlay-open");
        closeButton.focus();

        requestAnimationFrame(() => {
            renderLargeStructure(smiles, enlargedImage);
        });
    });

    closeButton.addEventListener("click", closeOverlay);

    overlay.addEventListener("click", event => {
        if (event.target === overlay) {
            closeOverlay();
        }
    });

    document.addEventListener("keydown", event => {
        if (event.key === "Escape" && !overlay.hidden) {
            closeOverlay();
        }
    });
}

interface SmiDrawerGlobal {
    apply(): void;
}

function applySmiDrawer(): void {
    const smiDrawer = (
        window as typeof window & {
            SmiDrawer?: SmiDrawerGlobal;
        }
    ).SmiDrawer;

    if (!smiDrawer) {
        console.error("SmiDrawer is not available on window");
        return;
    }

    const pendingImages = Array.from(
        document.querySelectorAll<HTMLImageElement>(
            "#ccs-table-body img[data-smiles]"
        )
    );

    if (pendingImages.length === 0) {
        return;
    }

    smiDrawer.apply();

    // SmiDrawer has already read the attributes, so these images will
    // not be selected and rendered again by the next apply() call.
    pendingImages.forEach(image => {
        image.removeAttribute("data-smiles");
    });
}

function renderImage(smiles: string): HTMLTableCellElement {
    const imageArea = document.createElement("td");
    imageArea.className = "preview-image-cell";

    const rowDetails = document.createElement("div");
    rowDetails.className = "ccs-row-details";

    const tableMolecule = document.createElement("div");
    tableMolecule.className = "ms1-table-molecule";

    const imageContainer = document.createElement("div");
    imageContainer.className = "im-container";

    const image = document.createElement("img");
    image.dataset.sourceSmiles = smiles;

    // setAttribute safely handles SMILES containing =, +, brackets, etc.
    image.setAttribute("data-smiles", smiles);
    image.setAttribute(
        "data-smiles-options",
        JSON.stringify({
            width: 150,
            height: 150,
            padding: 0.0
        })
    );
    image.setAttribute("data-smiles-theme", "dark");

    imageContainer.appendChild(image);
    tableMolecule.appendChild(imageContainer);
    rowDetails.appendChild(tableMolecule);
    imageArea.appendChild(rowDetails);

    return imageArea;
}


function toggleCcsColumn(show: boolean): void {
    const table = document.getElementById("ccs-prediction-table");

    table?.querySelectorAll<HTMLElement>(".preview-image-cell").forEach(cell => {
        cell.hidden = !show;
    });
}

export function check_selected(selectID: string) {
    const mySelect = document.getElementById(selectID) as HTMLSelectElement;
    const selectedValue = mySelect.value;
    return selectedValue
}

interface CcsPredictionResult {
    predictedCcs: number;
    mz: number;
}

async function ccs_prediction(): Promise<void> {
    const smilesElement = document.getElementById(
        "ccs-input"
    ) as HTMLInputElement;

    const smiles = smilesElement.value.trim();
    const adduct = check_selected("dropdown-adducts");

    await predictMultiple([
        {
            name: "Unnamed",
            smiles,
            adduct
        }
    ]);
}

async function open_ccs_file() {
    console.log("multiple ccs")
    let ccsfile = await get_csv();
    console.log("ccsfile:", ccsfile);
    
}

export async function get_csv() {
    console.log("hello")
    let result: string = await open({directory: false, multiple: false, 
        filters: [{name: 'Supported Files', extensions:['csv', 'tsv', 'txt']}]}) as string;
    data(result);
    glbl_result = result
    return result
}

async function data(result: string) {
    if (result && result.length > 0) {
        
        showPopup();
        let delimiter = ",";
        let file_contents: string[][] = [[]];

        file_contents = await read_input_csv(result as string, delimiter, dataset_transposed, deleted_ms1_rows);

        add_filename_to_top_of_popup(result)
        add_to_popup_table(file_contents)
        addPopupDropdownElements(file_contents[0] ?? []);
        //addPopupContinueElements()
        console.log("Executed data function:", result)
        glbl_result = result;
    } 
}

function add_to_popup_table(file_contents) {
    const tbody = document.getElementById("ccs-popup-tbody")
        tbody!.innerHTML = ""
        let upper_limit = 50
        if (file_contents.length < 50) {
            upper_limit = file_contents.length
        }
        for (let i = 0; i < upper_limit; i++) {
            tbody!.insertAdjacentHTML("beforeend", add_popup_table_row(file_contents[i], i));
            let trashbutton = document.getElementById(`ccs-popup-trash-button${i}`)
            trashbutton!.addEventListener("click", delete_popup_table_row);
        }

}

function delete_popup_table_row(event: Event) {
    const clickedButtonId = (event.target as HTMLElement).parentElement!.id;
    let row_idx = clickedButtonId.replace("ccs-popup-trash-button", "");

    const row = document.getElementById("ccs-popup-row" + row_idx as string);
    if (row) {
        row.remove();
    }
    deleted_ms1_rows.push(parseInt(row_idx));
}

function add_popup_table_row(items: string[], i: number) {
    let template = `
    <tr id="ccs-popup-row${i}">
        <td><button class="ccs-popup-trash-button" id="ccs-popup-trash-button${i}"><ion-icon name="trash-outline"></ion-icon></button></td>
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

function add_filename_to_top_of_popup(result){
    let resultString: string | null = Array.isArray(result) ? result.join(", ") : result;
    let filenameWindow = document.getElementById("ccs-popup-filename-window-name");
    filenameWindow!.textContent = resultString;
}

export function showPopup() {
    const overlay = document.getElementById('ccs-overlay');
    const popup = document.getElementById('ccs-popup');

    overlay!.style.display = 'block';
    popup!.style.display = 'flex';
    
    const closePopupButton = document.getElementById("ccs-popup-cancel-button");
    closePopupButton?.addEventListener("click", () => hidePopup())

    const overlayDiv = document.getElementById("ccs-overlay");
    overlayDiv?.addEventListener("click", () => hidePopup());

    eventListenerModule.addPopupTransposeListener();
    eventListenerModule.addPopupChangeDelimiterListener();
    eventListenerModule.addPopupContinueListener();

    (document.getElementById("ccs-popup-delimiter-input") as HTMLInputElement).value = ","
}

function hidePopup() {
    console.log("hiding")
    const overlay = document.getElementById('ccs-overlay');
    const popup = document.getElementById('ccs-popup');

    overlay!.style.display = 'none';
    popup!.style.display = 'none';

    eventListenerModule.removePopupTransposeListener();
    eventListenerModule.removePopupChangeDelimiterListener();
    eventListenerModule.removePopupContinueListener()
}

const popupTransposeButton = document.getElementById("ccs-popup-transpose-button");
async function popupTransposeListener(_event: Event) {
    dataset_transposed = !dataset_transposed;

    const delimiter = (
        document.getElementById(
            "ccs-popup-delimiter-input",
        ) as HTMLInputElement
    ).value;

    const fileContents = await read_input_csv(
        glbl_result as string,
        delimiter,
        dataset_transposed,
        deleted_ms1_rows,
    );

    add_to_popup_table(fileContents);
    addPopupDropdownElements(fileContents[0] ?? []);
}

function getPopupColumnMapping(): CcsColumnMapping {
    const mapping: Partial<CcsColumnMapping> = {};

    const selects =
        document.querySelectorAll<HTMLSelectElement>(
            ".ccs-popup-column-role",
        );

    selects.forEach((select) => {
        const role = select.value as CcsColumnRole;
        const columnIndex = Number(
            select.dataset.columnIndex,
        );

        if (!Number.isInteger(columnIndex) || role === "") {
            return;
        }

        switch (role) {
            case "name":
                mapping.name = columnIndex;
                break;

            case "smiles":
                mapping.smiles = columnIndex;
                break;

            case "adduct":
                mapping.adduct = columnIndex;
                break;

            case "m/z":
                mapping.mz = columnIndex;
                break;
        }
    });

    if (mapping.smiles === undefined) {
        throw new Error(
            "Select one CSV column as smiles before continuing",
        );
    }

    return mapping as CcsColumnMapping;
}

const popupChangeDelimiterButton = document.getElementById("ccs-popup-change-delimiter-button");
async function popupChangeDelimiterListener(_event: Event) {
    let delimiter = (document.getElementById("ccs-popup-delimiter-input") as HTMLInputElement).value as string;
    let file_contents: string[][] = [[]];
    file_contents = await read_input_csv(glbl_result as string, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    
    add_to_popup_table(file_contents);
}

const popupContinueButton = document.getElementById("ccs-popup-continue-button");
async function popupContinueListener(_event: Event) {
    const delimiter = (
        document.getElementById(
            "ccs-popup-delimiter-input",
        ) as HTMLInputElement
    ).value;

    let mapping: CcsColumnMapping;

    try {
        mapping = getPopupColumnMapping();
    } catch (error) {
        console.error(error);
        alert(String(error));
        return;
    }

    hidePopup();

    try {
        await read_csv(
            glbl_result,
            deleted_ms1_rows,
            delimiter,
            mapping,
        );
    } catch (error) {
        console.error(
            "Multiple CCS prediction failed:",
            error,
        );
    }
}

const eventListenerModule = {
    //Continue
    addPopupContinueListener() {
        popupContinueButton!.addEventListener("click", popupContinueListener);
    },
    removePopupContinueListener() {
        popupContinueButton?.removeEventListener("click", popupContinueListener);
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

async function read_input_csv(result: string, delimiter: string, transpose: boolean, deleted_rows: number[]) {
    console.log("awaiting...")
    const parsed_csv: string[][] = await invoke("read_ms1_csv", {
                                                        path: result, 
                                                        typeofcsv: "ccs",
                                                        n: 10000000,
                                                        delimiter: delimiter, 
                                                        transpose: transpose, 
                                                        deletedrows: deleted_rows
                                                    });
    console.log(parsed_csv)
    return parsed_csv
}

async function read_csv(
    result: string | string[] | null,
    deletedRows: number[],
    delimiter: string,
    mapping: CcsColumnMapping,
): Promise<void> {
    if (typeof result !== "string" || result.length === 0) {
        throw new Error("No CSV file was selected");
    }

    const tableBody =
        document.getElementById("ccs-table-body");

    if (!tableBody) {
        throw new Error("Could not find #ccs-table-body");
    }

    tableBody.innerHTML = "";

    const parsedCsv = await invoke<string[][]>(
        "read_ms1_csv",
        {
            path: result,
            typeofcsv: "ccs",
            n: 1_000_000,
            delimiter,
            transpose: !dataset_transposed,
            deletedrows: deletedRows,
        },
    );

    const smilesValues =
        parsedCsv[mapping.smiles] ?? [];

    const nameValues =
        mapping.name === undefined
            ? []
            : parsedCsv[mapping.name] ?? [];

    const adductValues =
        mapping.adduct === undefined
            ? []
            : parsedCsv[mapping.adduct] ?? [];

    const mzValues =
        mapping.mz === undefined
            ? []
            : parsedCsv[mapping.mz] ?? [];

    const molecules: CcsPredictionInput[] = [];

    for (
        let index = 0;
        index < smilesValues.length;
        index++
    ) {
        const smiles =
            String(smilesValues[index] ?? "").trim();

        // Skip an empty row or CSV header row.
        if (
            !smiles ||
            smiles.toLowerCase() === "smiles"
        ) {
            continue;
        }

        const name =
            String(nameValues[index] ?? "").trim() ||
            `Molecule ${index + 1}`;

        const adduct =
            String(adductValues[index] ?? "").trim() ||
            "[M+H]+";

        const sourceMzText =
            String(mzValues[index] ?? "").trim();

        const sourceMzNumber =
            sourceMzText === ""
                ? undefined
                : Number(sourceMzText);

        molecules.push({
            name,
            smiles,
            adduct,
            sourceMz:
                sourceMzNumber !== undefined &&
                Number.isFinite(sourceMzNumber)
                    ? sourceMzNumber
                    : undefined,
        });
    }

    if (molecules.length === 0) {
        throw new Error(
            "The selected file contains no valid SMILES",
        );
    }

    await predictMultiple(molecules);
}


function fill_table(item: [
    string,
    string,
    string,
    number,
    number
]): void {
    const tableBody = document.getElementById(
        "ccs-table-body"
    ) as HTMLTableSectionElement | null;

    if (!tableBody) {
        console.error("Could not find #ccs-table-body");
        return;
    }

    const row = document.createElement("tr");

    const columns = [
        { value: item[0], className: "name-column" },
        { value: item[1], className: "smiles-column" },
        { value: item[2], className: "adduct-column" },
        { value: item[3].toFixed(4), className: "mz-column" },
        { value: item[4].toFixed(2), className: "ccs-column" }
    ];

    columns.forEach(({ value, className }) => {
        const cell = document.createElement("td");
        cell.textContent = String(value);
        cell.className = className;
        row.appendChild(cell);
    });

    row.dataset.searchText = columns
        .map(column => String(column.value))
        .join(" ")
        .toLowerCase();

    const imageArea = renderImage(item[1]);
    imageArea.classList.add("structure-column");

    const previewSwitch = document.getElementById(
        "toggle-ccs-column"
    ) as HTMLInputElement | null;

    // A newly added preview inherits the current switch state.
    imageArea.hidden = previewSwitch
        ? !previewSwitch.checked
        : false;

    row.appendChild(imageArea);
    tableBody.appendChild(row);

    const currentSearch = (
        document.getElementById("ccs-table-search") as HTMLInputElement | null
    )?.value ?? "";

    row.style.display = rowMatchesCcsSearch(row, currentSearch)
        ? ""
        : "none";

    setTimeout(applySmiDrawer, 0);
}

interface SmiDrawerOptions {
    width?: number;
    height?: number;
    padding?: number;
}

interface SmiDrawerInstance {
    draw(
        smiles: string,
        target: HTMLImageElement,
        theme?: string
    ): void;
}

interface SmiDrawerGlobal {
    new(options?: SmiDrawerOptions): SmiDrawerInstance;
    apply(): void;
}

function renderLargeStructure(
    smiles: string,
    target: HTMLImageElement
): void {
    const SmiDrawerClass = (
        window as typeof window & {
            SmiDrawer?: SmiDrawerGlobal;
        }
    ).SmiDrawer;

    if (!SmiDrawerClass) {
        console.error("SmiDrawer is not available on window");
        return;
    }

    target.removeAttribute("src");

    const drawer = new SmiDrawerClass({
        width: 800,
        height: 800,
        padding: 0.05
    });

    drawer.draw(smiles, target, "dark");
}

/*
function removeItemsByIndices(arr: string[], indicesToRemove: number[]): string[] {
    return arr.filter((_, index) => !indicesToRemove.includes(index));
}


export function renderCCSHTML() {
    let ccs_div = document.getElementById("ccs-smiles-input");
}
*/

function guessColumnRole(columnName: string): CcsColumnRole {
    const normalized = columnName
        .trim()
        .toLowerCase()
        .replaceAll("_", " ");

    if (
        normalized === "name" ||
        normalized === "molecule name" ||
        normalized === "compound name"
    ) {
        return "name";
    }

    if (
        normalized === "smiles" ||
        normalized === "canonical smiles" ||
        normalized === "isomeric smiles"
    ) {
        return "smiles";
    }

    if (
        normalized === "adduct" ||
        normalized === "ion" ||
        normalized === "ion type"
    ) {
        return "adduct";
    }

    if (
        normalized === "m/z" ||
        normalized === "mz" ||
        normalized === "mass to charge"
    ) {
        return "m/z";
    }

    return "";
}

function addPopupDropdownElements(
    csvHeaderValues: string[],
): void {
    const table = document.getElementById(
        "ccs-popup-datatable",
    ) as HTMLTableElement | null;

    if (!table) {
        return;
    }

    const headerRow =
        table.tHead?.rows[0] ??
        table.rows[0];

    if (!headerRow) {
        return;
    }

    const roles: CcsColumnRole[] = [
        "",
        "name",
        "smiles",
        "adduct",
        "m/z",
    ];

    const usedRoles = new Set<CcsColumnRole>();

    // Cell 0 is the delete/trash column.
    for (
        let cellIndex = 1;
        cellIndex < headerRow.cells.length;
        cellIndex++
    ) {
        const cell = headerRow.cells[cellIndex];
        const csvColumnIndex = cellIndex - 1;

        const existingSelect =
            cell.querySelector<HTMLSelectElement>(
                ".ccs-popup-column-role",
            );

        const previousRole =
            existingSelect?.value as CcsColumnRole | undefined;

        const sourceHeader =
            csvHeaderValues[csvColumnIndex] ??
            `Column ${csvColumnIndex + 1}`;

        let selectedRole: CcsColumnRole = "";

        if (previousRole && !usedRoles.has(previousRole)) {
            selectedRole = previousRole;
        } else {
            const guessedRole = guessColumnRole(sourceHeader);

            if (
                guessedRole !== "" &&
                !usedRoles.has(guessedRole)
            ) {
                selectedRole = guessedRole;
            }
        }

        if (selectedRole !== "") {
            usedRoles.add(selectedRole);
        }

        const select = document.createElement("select");

        select.className = "ccs-popup-column-role";
        select.dataset.columnIndex =
            String(csvColumnIndex);

        for (const role of roles) {
            const option = document.createElement("option");

            option.value = role;
            option.textContent = role;

            select.appendChild(option);
        }

        select.value = selectedRole;
        select.title =
            `${sourceHeader} — Column ${csvColumnIndex + 1}`;

        /*
         * Each role may only be assigned to one column.
         * Selecting "smiles", for example, clears "smiles"
         * from any other column.
         */
        select.addEventListener("change", () => {
            const selectedRole =
                select.value as CcsColumnRole;

            if (selectedRole === "") {
                return;
            }

            document
                .querySelectorAll<HTMLSelectElement>(
                    ".ccs-popup-column-role",
                )
                .forEach((otherSelect) => {
                    if (
                        otherSelect !== select &&
                        otherSelect.value === selectedRole
                    ) {
                        otherSelect.value = "";
                    }
                });
        });

        cell.replaceChildren(select);
    }
}



interface CcsPredictionInput {
    name: string;
    smiles: string;
    adduct: string;
    sourceMz?: number;
}

interface CcsPredictionResult {
    predictedCcs: number;
    mz: number;
}

//interface CcsPredictionOutput extends CcsPredictionResult {
//    smiles: string | null;
//    adduct: string | null;
//}

interface CcsPredictionOutput extends CcsPredictionInput {
    predictedCcs: number | null;
    mz: number | null;
    error?: string;
}

async function predictMultiple(
    molecules: CcsPredictionInput[]
): Promise<CcsPredictionOutput[]> {
    if (ccsPredictionRunning) {
        console.warn("A CCS prediction is already running");
        return [];
    }

    const clearButton = document.getElementById(
        "clear-ccs-table-button"
    ) as HTMLButtonElement | null;

    if (clearButton) {
        clearButton.disabled = true;
    }

    const progressBar = document.getElementById(
        "ccs-prediction-progress-bar"
    ) as HTMLProgressElement | null;

    const progressText = document.getElementById(
        "ccs-progress-text"
    );

    const status = document.getElementById(
        "ccs-prediction-status"
    );

    const pauseButton = document.getElementById(
        "pause-ccs-prediction-button"
    ) as HTMLButtonElement | null;

    const cancelButton = document.getElementById(
        "cancel-ccs-prediction-button"
    ) as HTMLButtonElement | null;

    const tableBody =
        document.getElementById("ccs-table-body");

    if (!tableBody) {
        throw new Error("Could not find #ccs-table-body");
    }

    const total = molecules.length;
    let completed = 0;

    const results: CcsPredictionOutput[] = [];

    ccsPredictionRunning = true;
    ccsPredictionPaused = false;
    ccsPredictionCancelled = false;

    if (pauseButton) {
        pauseButton.disabled = false;
        updatePauseButton(false);
    }

    if (cancelButton) {
        cancelButton.disabled = false;
    }

    if (status) {
        status.textContent = "Running";
    }

    if (progressBar) {
        progressBar.max = Math.max(total, 1);
        progressBar.value = 0;
    }

    if (progressText) {
        progressText.textContent = `0 / ${total} (0%)`;
    }

    startCcsElapsedTimer();

    try {
        for (const molecule of molecules) {
            if (ccsPredictionCancelled) {
                break;
            }

            await waitWhileCcsPaused();

            if (ccsPredictionCancelled) {
                break;
            }

            let output: CcsPredictionOutput;

            try {
                const prediction =
                    await invoke<CcsPredictionResult>(
                        "sidecar_ccs",
                        {
                            smiles: molecule.smiles,
                            adduct: molecule.adduct
                        }
                    );

                output = {
                    ...molecule,
                    predictedCcs: prediction.predictedCcs,
                    mz: prediction.mz
                };
            } catch (error) {
                console.error(
                    `Prediction failed for ${molecule.smiles}:`,
                    error
                );

                output = {
                    ...molecule,
                    predictedCcs: null,
                    mz: null,
                    error: String(error)
                };
            }

            // Cancel may have been pressed while invoke() was running.
            if (ccsPredictionCancelled) {
                break;
            }

            results.push(output);

            if (
                output.mz !== null &&
                output.predictedCcs !== null
            ) {
                fill_table([
                    output.name,
                    output.smiles,
                    output.adduct,
                    output.mz,
                    output.predictedCcs
                ]);
            }

            completed += 1;

            const percentage =
                total === 0
                    ? 100
                    : Math.round((completed / total) * 100);

            if (progressBar) {
                progressBar.value = completed;
            }

            if (progressText) {
                progressText.textContent =
                    `${completed} / ${total} ` +
                    `(${percentage}%)`;
            }
        }

        const finalPercentage =
            total === 0
                ? 100
                : Math.round((completed / total) * 100);

        if (ccsPredictionCancelled) {
            if (progressText) {
                progressText.textContent =
                    `${completed} / ${total} ` +
                    `(${finalPercentage}%) — Cancelled`;
            }

            if (status) {
                status.textContent = "Cancelled";
            }
        } else if (status) {
            status.textContent = "Complete";
        }

        // Render molecule images once, after prediction stops.
        requestAnimationFrame(() => {
            applySmiDrawer();
        });

        return results;
    } finally {
        stopCcsElapsedTimer();

        ccsPredictionRunning = false;
        ccsPredictionPaused = false;
        ccsResumeResolver = null;

        if (pauseButton) {
            pauseButton.disabled = true;
            updatePauseButton(false);
        }

        if (cancelButton) {
            cancelButton.disabled = true;
        }
        if (clearButton) {
            clearButton.disabled = false;
        }
    }
}


export function createTableRowHTML(name: string, smiles: string, adduct: string, include_ccs: boolean): string {
  if (include_ccs == true) {
    return `
      <tr>
        <td class="clickable-cell">${name}</td>
        <td class="clickable-cell">${smiles}</td>
        <td class="clickable-cell">${adduct}</td>
        <td class="clickable-cell"></td>
        <td class="clickable-cell"></td>
        <td class="clickable-cell preview-column"></td>
          
      </tr>`;
  } else {
    return `
      <tr>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell preciew-column"></td>
      </tr>`;
  }
  
}

type CcsColumnRole = 
    | ""
    | "name"
    | "smiles"
    | "adduct"
    | "m/z";

interface CcsColumnMapping {
    smiles: number;
    name?: number;
    adduct?: number;
    mz?: number;
}
/*
function escapeCsvValue(value: string): string {
    // Quotes are doubled according to the CSV format.
    return `"${value.replaceAll('"', '""')}"`;
}
*/


interface CcsCsvRow {
    name: string;
    smiles: string;
    adduct: string;
    mz: string;
    predictedCcs: string;
}

async function saveCcsPredictions(): Promise<void> {
    const tableBody = document.getElementById("ccs-table-body");

    if (!tableBody) {
        console.error("CCS results table was not found");
        return;
    }

    const rows: CcsCsvRow[] = Array.from(
        tableBody.querySelectorAll<HTMLTableRowElement>("tr")
    ).map((row) => {
        const cells = row.querySelectorAll<HTMLTableCellElement>("td");

        return {
            name: cells[0]?.textContent?.trim() ?? "",
            smiles: cells[1]?.textContent?.trim() ?? "",
            adduct: cells[2]?.textContent?.trim() ?? "",
            mz: cells[3]?.textContent?.trim() ?? "",
            predictedCcs: cells[4]?.textContent?.trim() ?? "",
        };
    });

    if (rows.length === 0) {
        console.warn("There are no predictions to save");
        return;
    }

    try {
        const savedPath = await invoke<string | null>(
            "save_ccs_predictions",
            { rows }
        );

        if (savedPath) {
            console.log(`CCS predictions saved to: ${savedPath}`);
        }
    } catch (error) {
        console.error("Failed to save CCS predictions:", error);
    }
}

/*
document
    .getElementById("save-ccs-prediction-button")
    ?.addEventListener("click", () => {
        void saveCcsPredictions();
    });
*/