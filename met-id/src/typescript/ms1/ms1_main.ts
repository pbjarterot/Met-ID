import { invoke } from '@tauri-apps/api';
import { get_csv } from './ms1_popup';
import { identify, get_adjusted_ms1 } from './ms1_sidebar';
import { fill_dropdown, fill_under_dropdown } from "../dropdown";
import "./ms1_mass_error";
import "./ms1_table";
import "./ms1_popup";
import { addSearchbarListener } from "./ms1_searchbar";
import {add_matrix, add_metabolite, add_functional_group} from "./ms1_add_buttons";
import { convertTableToCSV } from './ms1_io';


window.addEventListener("DOMContentLoaded", async () => {
    document.getElementById("ms1-sidebar-add-metabolite")!.addEventListener("click", async () => add_metabolite());
    document.getElementById("ms1-sidebar-add-matrix")!.addEventListener("click", async () => add_matrix());
    document.getElementById("ms1-sidebar-add-functional-group")!.addEventListener("click", async () => add_functional_group());

    fill_dropdown(["HMDB (All)", "HMDB (Brain)", "HMDB (CSF)", "Lipidmaps"], "metabolome-dropdown");
    fill_dropdown(["Positive mode", "Negative mode", "FMP-10", "AMPP"], "matrix-dropdown");

    fill_under_dropdown.metabolites("metabolome-dropdown", "metabolome-checkbox-container")
    fill_under_dropdown.matrices("matrix-dropdown", "matrix-checkbox-container")
    

    document.getElementById("ms1-sidebar-open-file-button")!.addEventListener("click", async () => get_csv());
    document.getElementById("ms1-sidebar-identify-button2")!.addEventListener("click", async () => identify());
    document.getElementById("ms1-sidebar-export-adjusted")!.addEventListener("click", async () => get_adjusted_ms1());

    document.getElementById("matrix-dropdown")!.addEventListener("change", async () => check_options_for_sql_counter());
    document.getElementById("metabolome-dropdown")!.addEventListener("change", async () => check_options_for_sql_counter());

    var slider = document.getElementById("mzWindow") as HTMLInputElement;
    var output = document.getElementById("demo") as HTMLElement;
    output!.innerHTML = slider!.value + "ppm";

    slider!.oninput = function() {
        output!.innerHTML = (this as HTMLInputElement).value + " ppm";
    }


    document.getElementById("ms1-sidebar-export-button")!.addEventListener("click", () => {
        convertTableToCSV("ms1-datatable");
    })

    addSearchbarListener();
})





async function check_options_for_sql_counter() {
    fill_under_dropdown.metabolites("metabolome-dropdown", "metabolome-checkbox-container")
    fill_under_dropdown.matrices("matrix-dropdown", "matrix-checkbox-container")


    let met_selected: string = check_selected("metabolome-dropdown");
    let matrix_selected: string = check_selected("matrix-dropdown");

    let met_type: string[] = check_checkboxes("ms1-metabolome-div1");
    let adducts: string[] = check_checkboxes("ms1-metabolome-div2");

    let count: number = await invoke("sql_counter_tauri", {met: met_selected, mat: matrix_selected, typ: met_type, adducts:adducts});
    let db_size = document.getElementById("db_size");
    
    if (db_size !== null) {
        db_size.textContent = count.toString();
    }
}

export function check_checkboxes(divID: string) {
    const myDiv = document.getElementById(divID);
    const checkboxes = myDiv!.querySelectorAll("input[type='checkbox']");
    const checkedCheckboxes: HTMLInputElement[] = [];

    let checked_list:string[] = [];

    for (let i = 0; i < checkboxes.length; i++) {
        const checkbox = checkboxes[i] as HTMLInputElement;
        if (checkbox.checked) {
            checkedCheckboxes.push(checkbox);
            checked_list.push(checkbox.value)
        }
    }
    return checked_list
}

export function check_selected(selectID: string) {
    const mySelect = document.getElementById(selectID) as HTMLSelectElement;
    const selectedValue = mySelect.value;
    return selectedValue
}


