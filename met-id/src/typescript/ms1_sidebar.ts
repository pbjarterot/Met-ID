import { invoke } from '@tauri-apps/api';
import { open } from '@tauri-apps/api/dialog';
import { append_to_table } from './ms1_table';



window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("ms1-sidebar-open-file-button")!.addEventListener("click", async () => get_csv())
    document.getElementById("ms1-sidebar-identify-button2")!.addEventListener("click", async () => identify())
})

async function get_csv() {
    const result = await open({
        directory: false,
        multiple: false,
        filters: [{name: 'CSV Files', extensions:['csv']}],
    });
    
    if (result && result.length > 0) {
        const filePath = result;
        let parsed_csv: number[] = await invoke("parse_ms1_csv", {path: filePath});

        for (let i = 0; i < parsed_csv.length; i++) {
            await append_to_table(parsed_csv[i])
        }
    }
}

function check_checkboxes(divID: string) {
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
    //console.log(checkedCheckboxes);
    return checked_list
}

function check_selected(selectID: string) {
    const mySelect = document.getElementById(selectID) as HTMLSelectElement;
    const selectedValue = mySelect.value;
    return selectedValue
}




function get_ms1_input_peaks() {
    // Get a reference to the table element
    const table = document.getElementById("ms1-datatable") as HTMLTableElement;

    // Get a reference to the cells in the first column
    const cells = table.querySelectorAll("td:nth-child(2)");

    // Loop through the cells and extract their values
    const columnValues = [] as string[];
    for (let i = 0; i < cells.length; i++) {
        const myVariable: string | null = cells[i].textContent;
        if (myVariable !== null) {
            columnValues.push(myVariable);
        }
    }
    console.log(columnValues);
    return columnValues
}

function fill_ms1_results() {
    const table = document.getElementById("ms1-datatable") as HTMLTableElement;
    const tbody = table.querySelector('tbody');

    // Remove existing rows from the table body
    while (tbody!.firstChild) {
    tbody!.removeChild(tbody!.firstChild);
    }
    //let resultRow = ms1_rsults_data[0];
    // Loop through each row
    for (let i = 0; i < ms1_results_data.length; i++) {
        let oMass = "";
        let aMass = "";
        let names = "";
        let dMass = "";
        let dPPM = "";
        let tMass = "";
        let ms2 = "";
        let smiles_ = "";
        let smiles_for_img: string[] = [];
        let names_for_img: string[] = [];

        if (ms1_results_data[i].length === 1) {
            oMass = ms1_results_data[i][0].mz;
            aMass = ms1_results_data[i][0].mz;
            dMass = ms1_results_data[i][0].mz;
            tMass = ms1_results_data[i][0].mz;
            names = ms1_results_data[i][0].names;
            names_for_img.push(names);
            dPPM = "0";
            ms2 = "not-available"
            smiles_ = ms1_results_data[i][0].smiles;
            smiles_for_img.push(smiles_);
        } else {

            for (let j = 0; j < ms1_results_data[i].length; j++) {
                if (j < ms1_results_data[i].length) {
                    oMass += ms1_results_data[i][j].mz + "<br>";
                    aMass += ms1_results_data[i][j].mz + "<br>";
                    dMass += ms1_results_data[i][j].mz + "<br>";
                    tMass += ms1_results_data[i][j].mz + "<br>";
                    names += ms1_results_data[i][j].names + "<br>";
                    names_for_img.push(ms1_results_data[i][j].names);
                    dPPM += "0" + "<br>";
                    ms2 += `<span class="available"></span>` + "<br>";
                    smiles_ += ms1_results_data[i][j].smiles;
                    smiles_for_img.push(ms1_results_data[i][j].smiles);
                } else {
                    oMass += ms1_results_data[i][j].mz;
                    aMass += ms1_results_data[i][j].mz;
                    dMass += ms1_results_data[i][j].mz;
                    tMass += ms1_results_data[i][j].mz;
                    names += ms1_results_data[i][j].names;
                    names_for_img.push(ms1_results_data[i][j].names);
                    dPPM += "0";
                    ms2 = `<span class="available"></span>`;
                    smiles_ += ms1_results_data[i][j].smiles;
                    smiles_for_img.push(ms1_results_data[i][j].smiles);
                }
            }
        }


        const new_top_row = document.createElement("tr");
        new_top_row.setAttribute("class", "data");
        new_top_row.innerHTML = `
        <td><input type="checkbox" /></td>
        <td class="clickable-cell">${oMass}</td>
        <td class="clickable-cell">${aMass}</td>
        <td class="clickable-cell">${names}</td>
        <td class="clickable-cell">${dMass}</td>
        <td class="clickable-cell">${dPPM}</td>
        <td class="clickable-cell">${tMass}</td>
        <td class="clickable-cell">${ms2}</td>`
        tbody!.appendChild(new_top_row);


        const new_bottom_row = document.createElement("tr");
        new_bottom_row.setAttribute("class", "hidden-row");

        const image_area = document.createElement("td");
        image_area.setAttribute("colspan", "8");

        const row_details = document.createElement("div");
        row_details.setAttribute("class", "row-details");

        for (let index = 0; index < names_for_img.length; index++) {
            const table_molecule = document.createElement("div");
            table_molecule.setAttribute("class", "ms1-table-molecule");
            //table_molecule.classList.add("ms1-molecule")

            const imcontainer = document.createElement('div');
            imcontainer.setAttribute("class", "im-container");
            imcontainer.innerHTML = `<img data-smiles=${smiles_for_img[index]} data-smiles-options="{'width': 250, 'height': 250, 'padding': 0.0 }" data-smiles-theme='dark' />`

            const imcontainer2 = document.createElement('div');
            imcontainer2.setAttribute("class", "im-container2");
            
            imcontainer2.innerHTML += `<p class="ms1-molecule-name">${names_for_img[index]}</p>`

            const script = document.createElement('script');
            script.innerHTML = `
                        SmiDrawer.apply();
                        `;
            imcontainer.appendChild(script);
            table_molecule.append(imcontainer);
            table_molecule.append(imcontainer2);
            row_details.appendChild(table_molecule);
        }
        image_area.appendChild(row_details);

        new_bottom_row.appendChild(image_area);
        tbody!.appendChild(new_bottom_row);
    }    
    add_row_listeners();
}


function add_row_listeners() {
    const paragraphs = document.querySelectorAll('.im-container2 p');
    const heightThreshold = 10; // Adjust the height threshold as needed
  
    paragraphs.forEach((paragraph: Element) => {
      if ((paragraph as HTMLElement).offsetHeight > heightThreshold) {
        (paragraph as HTMLElement).classList.add('long-text');
      }
    });
    const clickableCells = document.querySelectorAll('.clickable-cell');

  clickableCells.forEach((cell) => {
    cell.addEventListener('click', () => {
      const row = cell.parentNode as HTMLTableRowElement;
      row.classList.toggle('clicked'); // Toggle the 'clicked' class on the parent row
    });
  });
}





























/*
    let j = 0;
    for (let i = 1; i < table.rows.length; i++) {
            if (i%2 === 1) {
            let matched_masses: string = "";
            let ppm_diff: string = "";
            let matched_name: string = "";

            console.log(i, j);
            for (let j=0; j < ms1_rsults_data[(i-1)/2].length; j++) {
                let ret_string: string = "";

                if (j > 0) {
                    ret_string += " <br> "
                }
                matched_masses += ret_string + ms1_rsults_data[(i-1)/2][j].mz
                //matched_masses += ret_string + ms1_rsults_data[i-1][j].mz
                matched_name += ret_string + ms1_rsults_data[(i-1)/2][j].names

            }
            const row = table.rows[i];

            // Input the column value by identifier
            const col1 = row.cells[2] as HTMLTableCellElement;
            const col2 = row.cells[3] as HTMLTableCellElement;
            const col3 = row.cells[4] as HTMLTableCellElement;

            // Update the column value
            col1.innerHTML = matched_masses;
            col2.innerHTML = matched_name;
            //col3.innerHTML = "Value 3";
        }
        else if (i%2===0) {
            const row = table.rows[i];
            const template = document.createElement("template");
            template.innerHTML =  `
            <td colspan="8">
              <div class="row-details">
                <div class="ms1-table-molecule">
                  <div class="im-container">
                    <img data-smiles=${"COOCCCCCCO"} data-smiles-options= "{'width': 200, 'height': 160, 'padding': 0.0 }"  data-smiles-theme='dark' />
                  </div>
                  <div class="im-container2">
                    <p class="ms1-molecule-name">${"Something"}</p>
                  </div>
                  
                </div>
                <div class="ms1-table-molecule">
                  <div class="im-container">
                    <img data-smiles="CCOCCCCC" data-smiles-options="{'width': 200, 'height': 160, 'padding': 0.0 }"  data-smiles-theme='dark' />
                  </div>
                  <div class="im-container2">
                    <p class="ms1-molecule-name">Dopamine</p>
                  </div>
                </div>
              </div>
            </td>`
            const script = document.createElement('script');
            script.innerHTML = `
                        SmiDrawer.apply();
                        `;
            const templateContent = template.content;
            templateContent.appendChild(script);

            row.innerHTML = "";
            row.appendChild(templateContent);
        }

    }
}
*/



let ms1_results_data: Array<Array<Record<string, string>>>;

async function identify() {
    let met_selected: string = check_selected("metabolome-dropdown");
    let matrix_selected: string = check_selected("matrix-dropdown");

    let met_type: string[] = check_checkboxes("ms1-metabolome-div1");
    let adducts: string[] = check_checkboxes("ms1-metabolome-div2");

    const dictionary = {
    "Carboxylic": "caboxylicacids",
    "Primary": "primaryamines",
    }
    console.log(adducts);

    for (let i = 0; i < adducts.length; i++) {
        adducts[i] = dictionary[adducts[i]]?.value ?? adducts[i]
    }

    let mass_error_input: string = "0.0000";//(document.getElementById("ms1-error-input") as HTMLInputElement)!.value as string;
    let input_masses: string[] = get_ms1_input_peaks();


    ms1_results_data = await invoke("sql_handler", {met: met_selected, mat: matrix_selected, typ: met_type, 
                                         adducts:adducts, massError: mass_error_input, masses: input_masses});

    fill_ms1_results();
}
