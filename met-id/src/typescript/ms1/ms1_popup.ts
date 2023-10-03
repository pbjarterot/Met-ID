import { fill_dropdown } from '../dropdown';
import { invoke } from '@tauri-apps/api';
import { open } from '@tauri-apps/api/dialog';
import { createTableRowHTML } from './ms1_table';

let deleted_ms1_rows: number[] = [];
let dataset_transposed: boolean = false;
let glbl_result: string | string[] | null = "";
let inputType: string = "";


async function data(result: string, result_type: string) {
    if (result && result.length > 0) {
        if (deleted_ms1_rows.length > 0) {
            deleted_ms1_rows.length = 0
        }
        
        showPopup();
        let delimiter = ";";
        let file_contents: string[][] = [[]];
        if (result_type === "filename") {
            inputType = "file";
            file_contents = await read_input_csv(result as string, delimiter, dataset_transposed, deleted_ms1_rows );
        } else {
            if (containsDelimiter(result)) {
                inputType = "ctrlv";
                file_contents = await read_ctrl_v(result as string, delimiter, dataset_transposed, deleted_ms1_rows )
            } else {
                return
            }
            
        }   
        console.log(file_contents);
        add_filename_to_top_of_popup(result)
        add_to_popup_table(file_contents)
        addPopupDropdownElements();
        //addPopupContinueElements()
        glbl_result = result;
    } 
}

export function get_ctrl_v_data(result:string) {
    data(result, "content");
}

export async function get_csv() {
    let result: string = await open({directory: false, multiple: false, filters: [{name: 'Delimiter Separated Files', extensions:['csv', 'tsv']}, {name: 'Text Files', extensions:['txt']}]}) as string;
    data(result, "filename");
    
}

function containsDelimiter(str: string): boolean {
    return str.indexOf(',') !== -1 || str.indexOf(';') !== -1 || str.indexOf('\t') !== -1;
}

async function read_input_csv(result: string, delimiter: string, transpose: boolean, deleted_rows: number[]) {
    const parsed_csv: string[][] = await invoke("read_ms1_csv", {path: result, delimiter: delimiter, transpose: transpose, deletedrows: deleted_rows});
    return parsed_csv
}

async function read_ctrl_v(result: string, delimiter: string, transpose: boolean, deleted_rows: number[]) {
    const parsed_csv: string[][] = await invoke("read_ms1_ctrlv", {content: result, delimiter: delimiter, transpose: transpose, deletedrows: deleted_rows});
    return parsed_csv
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

function addPopupDropdownElements() {
    let table = document.getElementById("ms1-popup-datatable") as HTMLTableElement;
    if (table) {
        // Access the first row of the table
        const headerRow = table.rows[0];

        // Array to store the column names
        const columnNames: string[] = [];

        // Iterate over the cells of the header row
        for (let i = 0; i < headerRow.cells.length; i++) {
            const thElement = headerRow.cells[i];

            // Extract the text content of the th element
            const columnName = thElement.textContent;

            // Add the column name to the array
            columnNames.push(columnName!);
        }
        fill_dropdown(columnNames, "ms1-popup-mz-column-dropdown")
        //fill_dropdown(columnNames, "ms1-popup-ccs-column-dropdown")

        const dropdown = document.getElementById('ms1-popup-mz-column-dropdown') as HTMLSelectElement;

        // Set the value of the dropdown to select a specific item
        const selectedItemValue = 'Column 1';
        dropdown.value = selectedItemValue;

        const inputElement = document.getElementById('ms1-popup-delimiter-input') as HTMLInputElement;

        // Set the value of the input element
        inputElement.value = ';';
    }

}

const popupTransposeButton = document.getElementById("ms1-popup-transpose-button");
async function popupTransposeListener(_event: Event) {
    dataset_transposed = !dataset_transposed;
    let delimiter = (document.getElementById("ms1-popup-delimiter-input") as HTMLInputElement).value as string;
    let file_contents: string[][] = [[]];
    if (inputType === "file") {
        file_contents = await read_input_csv(glbl_result as string, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
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
        file_contents = await read_input_csv(glbl_result as string, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    } else if (inputType === "ctrlv"){
        file_contents = await read_ctrl_v(glbl_result as string, delimiter, dataset_transposed, deleted_ms1_rows as number[]);
    }
    
    add_to_popup_table(file_contents);
}

const popupContinueButton = document.getElementById("ms1-popup-continue-button");
async function popupContinueListener(_event: Event) {
    let delimiter = (document.getElementById("ms1-popup-delimiter-input") as HTMLInputElement).value as string;
    read_csv(glbl_result, deleted_ms1_rows, delimiter);
    hidePopup();
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

async function read_csv(result: string| string[] | null, deleted_rows: number[], delimiter: string) {
    console.log("read_csv was used");
    const table_body = document.getElementById("ms1-table-body");
    table_body!.innerHTML = "";
    const selectedColumn = parseInt(((document.getElementById("ms1-popup-mz-column-dropdown") as HTMLSelectElement).value as unknown as string).replace("Column ", "")) -1 ;
    let parsed_csv: number[][] = [[]];
    if (inputType === "file") {
        parsed_csv = await invoke("read_ms1_csv", {path: result, delimiter: delimiter, transpose: !dataset_transposed, deletedrows: deleted_rows});
    } else if (inputType === "ctrlv"){
        parsed_csv = await invoke("read_ms1_ctrlv", {content: result, delimiter: delimiter, transpose: !dataset_transposed, deletedrows: deleted_rows});
    }
    //let parsed_csv: number[][] = await invoke("read_ms1_csv", {path: result, delimiter: delimiter, transpose: !dataset_transposed, deletedrows: deleted_rows});
    let parsed_csv2 = parsed_csv[selectedColumn]

    
    const parsed_csv3: number[] = removeItemsByIndices(parsed_csv2, deleted_rows)
    console.log(parsed_csv3);


    let allRowsHTML = '';

    for (let i = 0; i < parsed_csv3.length; i++) {
        allRowsHTML += createTableRowHTML(parsed_csv3[i]);
    }
    //console.log(allRowsHTML);
    // Append all rows in one go
    table_body!.insertAdjacentHTML('beforeend', allRowsHTML);
}
