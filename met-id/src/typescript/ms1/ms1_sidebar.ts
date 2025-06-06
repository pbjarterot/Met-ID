import { invoke } from '@tauri-apps/api/core';
import { check_checkboxes, check_selected } from './ms1_main';
import { createBottomRow } from './ms1_table';
import { saveCsv } from './ms1_io';

let imageData = {"names": [[]] as string[][], "smiles": [[]] as string[][], "hmdb": [[]] as string[][]}


export async function identify() {
    let id_res_string = document.getElementById("identification-results-string") as HTMLParagraphElement;
    id_res_string.textContent = "";
    let met_selected: string = check_selected("metabolome-dropdown");
    let matrix_selected: string = check_selected("matrix-dropdown");

    let met_type: string[] = check_checkboxes("ms1-metabolome-div1");
    let adducts: string[] = check_checkboxes("ms1-metabolome-div2");

    let mass_error_input: string = (document.getElementById("ms1-error-input-text") as HTMLInputElement)!.value as string;
    let massWindow: string = (document.getElementById("mzWindow") as HTMLInputElement)!.value as string;
    let input_masses: string[] = get_ms1_input_peaks();

    let append: boolean = (document.getElementById("ms1-append-results") as HTMLInputElement)!.checked;
    console.log(append);
    let startTime = performance.now();

    let ms1_results_data: Array<Array<Record<string, string>>> = await invoke("sql_handler_tauri", {met: met_selected, mat: matrix_selected, typ: met_type, adducts:adducts, massError: mass_error_input, masses: input_masses, mzwindow: massWindow});

    let endTime = performance.now();
    let timeTaken = endTime - startTime;
    console.log(`Time taken: ${timeTaken} milliseconds`)
    fill_ms1_results(ms1_results_data, append);

    count_identified_percent();
}

function count_identified_percent() {
    const table = document.getElementById("ms1-datatable") as HTMLTableElement;
    const rows = table.querySelectorAll('tr');
    var numRows: number = 0;
    var numIdentified: number = 0;
    rows.forEach((row) => {
        if (row.className === "data") {
            numRows += 1;
            const cells = row.querySelectorAll('td, th');
            if (cells.length > 0 && cells[3].textContent != "" ) {
                numIdentified += 1
            }
        }
        
    });
    let id_res_string = document.getElementById("identification-results-string") as HTMLParagraphElement;
    id_res_string.textContent = `Out of the ${String(numRows)} m/z features,\n ${numIdentified} (${Math.round((numIdentified/numRows)*1000)/10}%) have been identfied`;
    console.log(numRows, numIdentified, numIdentified/numRows)
}

function parse_cell_text(cell: Element) {
    var innerHTML = cell!.innerHTML;

    // Find the index of the first occurrence of "<br>"
    var brIndex = innerHTML.indexOf("<br>");
    

    // Find the index of the first occurrence of "<br>"
    var brIndex = innerHTML.indexOf("<br>");

    if (brIndex !== -1) {
        // Remove everything starting from the first <br> tag
        var updatedInnerHTML = innerHTML.substring(0, brIndex);

        // Trim any leading or trailing whitespace
        updatedInnerHTML = updatedInnerHTML.trim();

        // Display the updated inner HTML
        return updatedInnerHTML;
    } else {
        return innerHTML;
    }
}

export async function get_adjusted_ms1() {
    let inputMasses: string[] = get_ms1_input_peaks();
    let massErrorInput: string = (document.getElementById("ms1-error-input-text") as HTMLInputElement)!.value as string;
    let adjusted_ms1s: Record<string, string>[] = await invoke("calculate_adjusted_mass", {masses:inputMasses, massError: massErrorInput});
    const dataRows = adjustedProcessDataGroup(adjusted_ms1s);
    saveCsv(dataRows.oMass, dataRows.aMass);
}

function get_ms1_input_peaks() {
    // Get a reference to the table element
    const table = document.getElementById("ms1-datatable") as HTMLTableElement;

    // Get a reference to the cells in the first column
    const cells = table.querySelectorAll("td:nth-child(2)");

    // Loop through the cells and extract their values
    const columnValues = [] as string[];
    for (let i = 0; i < cells.length; i++) {
        const myVariable: string = parse_cell_text(cells[i]);
        if (myVariable !== null) {
            columnValues.push(myVariable);
        }
    }
    return columnValues
}

function scientific_format(input_string: string) {
    const number = parseFloat(input_string);
    const scientificNotation = number.toExponential(5);
    if (isNaN(number)) {
        return "";
    } else {
        return scientificNotation.toString();
    }
}

function formNumber(numberString: string): string {
    const number = parseFloat(numberString);
    const formattedNumber = number.toFixed(6);
    if (isNaN(number)) {
        return "";
    } else {
        return formattedNumber.toString();
    }
    
}

function removeExistingRows(tbody: HTMLTableSectionElement) {
    while (tbody.firstChild) {
        tbody.removeChild(tbody.firstChild);
    }
}

function createDataCells(data: string) {
    const td = document.createElement("td");
    td.className = "clickable-cell";
    td.innerHTML = data;
    return td;
}

async function renderImagesAsync(idx: number): Promise<HTMLElement> {
    const image_area = document.createElement("td");
    image_area.setAttribute("colspan", "11");
    console.log(imageData[idx])

    const row_details = document.createElement("div");
    row_details.setAttribute("class", "row-details");
    
    const renderPromises = imageData.names[idx].map((name, index) => new Promise<void>((resolve) => {
        if (imageData.smiles[idx][index].length === 0) {
            resolve();
            return;
        }

        const table_molecule = document.createElement("div");
        table_molecule.className = "ms1-table-molecule";

        const imcontainer = document.createElement("div");
        imcontainer.className = "im-container";

        imcontainer.innerHTML = `<img data-smiles=${imageData.smiles[idx][index]} data-smiles-options="{'width': 250, 'height': 250, 'padding': 0.0 }" data-smiles-theme='dark' />`;

        const imcontainer2 = document.createElement("div");
        imcontainer2.className = "im-container2";
        imcontainer2.innerHTML = `<p class="ms1-molecule-name" onclick="window.open('https://hmdb.ca/metabolites/${imageData.hmdb[idx][index]}', '_blank')">${name}</p> `; //<button onclick="window.open('https://hmdb.ca/metabolites/${hmdb_id[index]}', '_blank')">HMDB</button>`;

        setTimeout(() => {
            const script = document.createElement("script");
            script.innerHTML = `SmiDrawer.apply();`;
            imcontainer.appendChild(script);
            table_molecule.appendChild(imcontainer);
            table_molecule.appendChild(imcontainer2);
            row_details.appendChild(table_molecule);
            resolve();
        }, 0);
    }));

    await Promise.all(renderPromises);
    
    image_area.appendChild(row_details);
    return image_area;
}

async function fill_ms1_results(ms1_results_data: Array<Array<Record<string, string>>>, append: boolean) {
    //const table = document.getElementById("ms1-datatable") as HTMLTableElement;
    const tbody = document.getElementById("ms1-table-body") as HTMLTableSectionElement;

    // Remove existing rows from the table body
    if (!append) {
        removeExistingRows(tbody);
    }
    

    // Create a document fragment to batch DOM manipulations
    const fragment = document.createDocumentFragment();

    // Loop through the data and create rows
    ms1_results_data.forEach(async (resultGroup, index) => {
        const dataRows = processDataGroup(resultGroup);
    
            if (!append) {
                const topRow = createTopRow(dataRows, index);
                fragment.appendChild(topRow);
                imageData.names[index] = dataRows.names;
                imageData.smiles[index] = dataRows.smiles;
                imageData.hmdb[index] = dataRows.hmdb;

            topRow.addEventListener('click', async () => {
                topRow.classList.toggle('clicked'); // Toggle the 'clicked' class on the parent row
                let bottomRow = document.getElementById("hidden-row-" + index);
            
                // Check if topRow has the 'clicked' class
                if (topRow.classList.contains('clicked')) {
                    // When topRow is clicked for the first time
                    let image_area = await renderImagesAsync(index);
                    bottomRow!.appendChild(image_area);
                } else {
                    // When topRow is clicked again
                    bottomRow!.innerHTML = ''; // Clear the contents of bottomRow
                }
            });
            
            // Create the bottom row with images if names are present
            if (dataRows.names.length > 0) {
                const bottomRow = createBottomRow(index, dataRows.names, dataRows.smiles);
                fragment.appendChild(bottomRow);
                
            }
            
        }
        else {
            AppendTopRow(dataRows, index);
            console.log("append", imageData.names[index], dataRows.names)
            imageData.names[index].push(...dataRows.names);
            console.log("append", imageData.names[index], dataRows.names)
            imageData.smiles[index].push(...dataRows.smiles);
            imageData.hmdb[index].push(...dataRows.hmdb);
        }
    });

    // Append the fragment to the table body
    tbody.appendChild(fragment);

}

function AppendTopRow(dataRows: { checkbox: any; oMass: any; aMass: any; names: any; adduct: any; dMass: any; tMass: any; dPPM: any; ms2: any; smiles: string[]; names_for_img?: string[]; formula: any; coverage: any; hmdb: any}, index: number) {
    const dataCells = [
        dataRows.checkbox,
        dataRows.oMass,
        dataRows.aMass,
        dataRows.names,
        dataRows.adduct,
        dataRows.formula,
        dataRows.dMass,
        dataRows.dPPM,
        dataRows.tMass,
        dataRows.ms2,
        dataRows.coverage
    ];
    // Get the table row by ID
    const row = document.getElementById("data-" + index) as HTMLTableRowElement;
    if (row.cells[3].innerHTML === "") {
        return;
    }
    if (row.cells[3].innerHTML.includes(dataCells[3][0])) {
        return;
    }
    console.log()
    // Check if the row exists
    if (!row) {
        console.error('Row not found');
        return;
    }

    // Append new data to each cell
    dataCells.forEach((cellData, index) => {
        // Get the cell at the current index, or create it if it doesn't exist
        const cell = row.cells[index];
        console.log(cell);
        //cell.innerHTML += '<br>';
        // Append each piece of new data in the cell, separated by <br>
        cellData.forEach((data: string, _idx: number) => {
            cell.innerHTML += '<br>';
            cell.innerHTML += data;
        });
    });
}




function processDataGroup(resultGroup: Record<string, string>[]) {
    let checkbox: string[] = []
    let oMass: string[] = [];
    let aMass: string[] = [];
    let names: string[] = [];
    let adduct: string[] = [];
    let dMass: string[] = [];
    let tMass: string[] = [];
    let dPPM: string[] = [];
    let ms2: string[] = [];
    let smiles: string[] = [];
    let formula: string[] = [];
    let hmdb: string[] = []
    let coverage: string[] = []

    for (const result of resultGroup) {
        checkbox.push('<input type="checkbox" /></td>');
        oMass.push(formNumber(result.oMass));
        aMass.push(formNumber(result.aMass));
        dMass.push(scientific_format(result.dMass));
        tMass.push(formNumber(result.tMass));
        names.push(result.names);
        adduct.push(result.matrix);
        dPPM.push(formNumber(result.dPPM));
        ms2.push(`<span class="${result.msms}"></span>`);
        smiles.push(result.smiles);
        formula.push(result.formula);
        hmdb.push(result.accession);
        coverage.push(`${result.coverage} / ${result.possible_derivs}`);
    }

    return {
        checkbox,
        oMass,
        aMass,
        names,
        adduct,
        dMass,
        tMass,
        dPPM,
        ms2,
        smiles,
        formula,
        hmdb,
        coverage
    };
}

function adjustedProcessDataGroup(resultGroup: Record<string, string>[]) {
    let oMass: string[] = [];
    let aMass: string[] = [];

    for (const result of resultGroup) {
        oMass.push(formNumber(result.oMass));
        aMass.push(formNumber(result.aMass));
    }

    return {
        oMass,
        aMass,
    };
}

function createTopRow(dataRows: { checkbox: any; oMass: any; aMass: any; names: any; adduct: any; dMass: any; tMass: any; dPPM: any; ms2: any; smiles?: string[]; names_for_img?: string[]; formula: any; coverage: any}, index: number) {
    const row = document.createElement("tr");
    row.className = "data";
    row.id = "data-" + index;

    const dataCells = [
        dataRows.checkbox,
        dataRows.oMass,
        dataRows.aMass,
        dataRows.names,
        dataRows.adduct,
        dataRows.formula,
        dataRows.dMass,
        dataRows.dPPM,
        dataRows.tMass,
        dataRows.ms2,
        dataRows.coverage
    ];
    for (const cellData of dataCells) {
        row.appendChild(createDataCells(cellData.join("<br>")));
    }
    //console.log(row);
    return row;
    
}
