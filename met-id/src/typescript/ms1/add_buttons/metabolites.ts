import { invoke } from "@tauri-apps/api/core";
import { change_slidein, slideInDiv, generateCheckboxes, addToHTML, createAndAttach, 
         remove_warning_message, add_to_db_rust, add_warning_message, removeFromDBInstance, cancellationToken } from "../ms1_add_buttons";


let MetaboliteListener: ((e: MouseEvent) => void) | null = null;
let MetaboliteRefreshListener: ((e: MouseEvent) => void) | null = null;

let SubmitMetaboliteElement: HTMLButtonElement | null = null;
let MetaboliteRefreshElement: HTMLSpanElement | null = null;


export async function add_metabolite() {
  change_slidein();
	slideInDiv!.innerHTML = "";
	slideInDiv!.innerHTML = await addToHTML("add-metabolite-to-db", "Metabolite", [["Name", "add-name-to-db"], ["SMILES", "add-smiles-to-db"]]);

	if (slideInDiv?.classList.contains('slide-in--active')) {
		SubmitMetaboliteElement = document.getElementById("add-metabolite-to-db") as HTMLButtonElement;
		MetaboliteRefreshElement = document.getElementById("ms1-slidein-table-h1-refresh") as HTMLSpanElement;
		createAndAttach("metabolite");
		update_user_metabolites();
	}
}

export async function addToHTML_Metabolite(table:string) {
    let checkboxes = generateCheckboxes(["Endogenous", "Exogenous", "Unspecified"], "Type of Metabolite");

		let checkbox_labels: string[] = await invoke("get_tissues_tauri", {});
		checkboxes += generateCheckboxes(checkbox_labels, "Found in tissue:");

		table += `
			<div class="ms1-slidein-spreadsheet-div" id="ms1-slidein-spreadsheet-div">
			<h1> User metabolites 
				<span class="ms1-slidein-table-h1-refresh-icon" id="ms1-slidein-table-h1-refresh">
					<ion-icon name="refresh"></ion-icon>
				</span>
			</h1>
				<div class="ms1-slidein-spreadsheet-table">
				<table class="ms1-slidein-datatable" id="ms1-slidein-datatable">
					<thead>
						<tr>
							<th>ID</th>
							<th>Name</th>
							<th>Smiles</th>
							<th>Formula</th>
							<th>m/z</th>
							<th>  </th>
							<!-- Add more column headers as needed -->
						</tr>
					</thead>
					<tbody id="ms1-slidein-tbody">
						<!-- Add more rows as needed -->
					</tbody>
				</table>
				</div>
			</div>
			`;
	return [checkboxes, table];
}

interface USER_METABOLITES {
	ids: string[],
	names: string[],
	smiless: string[],
	formulas: string[],
	mzs: string[]
}

export async function update_user_metabolites() {
	let metabolites: USER_METABOLITES = await invoke("update_user_metabolites_tauri");

	await renderSmallItemsFor_UserMetabolites(metabolites[0], 
																						metabolites[1], 
																						metabolites[2], 
																						metabolites[3], 
																						metabolites[4])
}

export async function remove_row_from_user_metabolites(row_id: string) {
	const parts = row_id.split("-");
	const lastpart = parts[parts.length - 1];
	const rowid = parseInt(lastpart);
	const done: number = await invoke("remove_row_from_user_metabolites_tauri", {rowid:rowid});
	return done;
}

async function renderSmallItemsFor_UserMetabolites(ids: string[], names: string[], smiless: string[], formulas: string[], mzs: string[]) {
	const res_div = document.getElementById("ms1-slidein-tbody") as HTMLTableElement;
	res_div.innerHTML = "";
	let i = 0;
    let rowIds: string[] = []


	function appendBatch() {
		if (cancellationToken.isCancelled) {
			console.log("Operation cancelled");
			return;
		}
		const template = document.createElement("template");
		for (let j = 0; j < ids.length; j++) {
			//checking for 0eV is not optimal when searching, as fragments could only appear in e.g. 30eV
			if (cancellationToken.isCancelled) {
				console.log("Operation cancelled");
				return;

			}
			if (ids.length > 0) {
				template.innerHTML = `<tr>
										<td>${ids[i]}</td>
										<td>${names[i]}</td>
										<td>${smiless[i]}</td>
										<td>${formulas[i]}</td>
										<td>${mzs[i]}</td>
										<td>
												<button class="ms1-slidein-table-trashbutton" id="ms1-slidein-table-trashbutton-${ids[i]}">
														<ion-icon name="trash-outline"></ion-icon>
												</button>
										</td>
									</tr>`;
				const templateContent = template.content;
				res_div.appendChild(templateContent);
                rowIds.push(`ms1-slidein-table-trashbutton-${ids[i]}`)
			}
        i++
		}

        removeFromDBInstance.addListeners(rowIds, "metabolites")
		
		if (i < names.length-1) {
			setTimeout(appendBatch, 0);
		}
	}
	appendBatch();
}

export function createAndAttachMetabolites() {
  MetaboliteListener = async () => {
    let name = document.getElementById("add-name-to-db") as HTMLInputElement;
    let smiles = document.getElementById("add-smiles-to-db") as HTMLInputElement;

    const smilesResult = true;

    const TissueCheckboxObject = {};
    const MetTypeCheckboxObject = {};

    const divElement = document.getElementById('add-to-database-checkboxes');
    if (divElement) {
      const checkboxes = divElement.querySelectorAll('input[type=checkbox]');
      checkboxes.forEach((checkbox) => {
        const inputCheckbox = checkbox as HTMLInputElement;

        const endo_exo_ids = ["add-Endogenous", "add-Exogenous", "add-Unspecified"];
        if (endo_exo_ids.includes(inputCheckbox.id)) {
          MetTypeCheckboxObject[inputCheckbox.id.replace("add-", "")] = inputCheckbox.checked
        } else {
          TissueCheckboxObject[inputCheckbox.id.replace("add-", "")] = inputCheckbox.checked
        }
      });
    }
  
    if (smilesResult) {
      remove_warning_message();
      add_to_db_rust(name!.value, smiles!.value, "metabolite", MetTypeCheckboxObject, TissueCheckboxObject, [])
    } else {
      add_warning_message("SMILES");
    }
    update_user_metabolites();

  };
  SubmitMetaboliteElement?.addEventListener("click", MetaboliteListener);

  MetaboliteRefreshListener = async () => {
    update_user_metabolites()
  }
  MetaboliteRefreshElement?.addEventListener("click", MetaboliteRefreshListener);

}

export function destroyAndDetachMetabolite() {
  if (MetaboliteListener) {
		SubmitMetaboliteElement!.removeEventListener("click", MetaboliteListener);
		MetaboliteListener = null;
	}
	if (MetaboliteRefreshListener) {
		MetaboliteRefreshElement!.removeEventListener("click", MetaboliteRefreshListener);
		MetaboliteRefreshListener = null;
	}
}


