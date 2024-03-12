import { invoke } from "@tauri-apps/api";
import { addToHTML, add_to_db_rust, add_warning_message, cancellationToken, change_slidein, createAndAttach, generateCheckboxes, removeFromDBInstance, remove_warning_message, slideInDiv } from "../ms1_add_buttons";


let MatrixListener: ((e: MouseEvent) => void) | null = null;
let MatrixRefreshListener: ((e: MouseEvent) => void) | null = null;

let SubmitMatrixElement: HTMLButtonElement | null = null;
let MatrixRefreshElement: HTMLSpanElement | null = null;


function isNumeric(str: string): boolean {
	return !isNaN(Number(str));
}




export function destroyAndDetachMatrix() {
  if (MatrixListener) {
		SubmitMatrixElement!.removeEventListener("click", MatrixListener);
		MatrixListener = null;
	}
	if (MatrixRefreshListener) {
		MatrixRefreshElement!.removeEventListener("click", MatrixRefreshListener);
		MatrixRefreshListener = null;
	}
}




export function createAndAttachMatrices() {
  console.log("Hello world");
		MatrixListener = async () => {
			let name = document.getElementById("add-name-to-db") as HTMLInputElement;
			let smiles_mz = document.getElementById("add-smiles-mz-to-db") as HTMLInputElement;
			let smilesResult = true;//await check_smiles(smiles_mz!.value);
			console.log(smilesResult);

			const matrixCheckboxObject = {}

			const divElement = document.getElementById('add-to-database-checkboxes');
			if (divElement) {
				const checkboxes = divElement.querySelectorAll('input[type=checkbox]');
				checkboxes.forEach((checkbox) => {
					const inputCheckbox = checkbox as HTMLInputElement;
					matrixCheckboxObject[inputCheckbox.id.replace("add-", "")] = inputCheckbox.checked
				});
			}

			const inputs = document.querySelectorAll('.add-to-database-adducts input');
			const values: string[] = [];

			inputs.forEach((input) => {
				values.push((<HTMLInputElement>input).value);
			});


			const matrixAdductsObject = {}

			for (let i = 0; i < inputs.length; i++) {
				let i_obj = inputs[i].id
				matrixAdductsObject[i_obj] = values[i];
			}


			if (smilesResult && !isNumeric(smiles_mz!.value)) {
				remove_warning_message();
				add_to_db_rust(name!.value, smiles_mz!.value, "matrix", matrixCheckboxObject, {}, matrixAdductsObject)
				console.log(name!.value, smiles_mz!.value, "matrix", matrixCheckboxObject, {}, matrixAdductsObject)
			} else if (isNumeric(smiles_mz!.value)) {
				remove_warning_message();
				add_to_db_rust(name!.value, smiles_mz!.value, "matrix", matrixCheckboxObject, {}, matrixAdductsObject)
			} else {
				add_warning_message("Value is not a number and SMILES")
			}


		};
		SubmitMatrixElement!.addEventListener("click", MatrixListener);
		MatrixRefreshListener = async () => {
			update_user_matrices()
		}
		MatrixRefreshElement?.addEventListener("click", MatrixRefreshListener);

}



export async function add_matrix() {
  change_slidein();
  slideInDiv!.innerHTML = "";
slideInDiv!.innerHTML = await addToHTML("add-matrix-to-db", "Matrix",[["Name", "add-name-to-db"], ["SMILES or m/z", "add-smiles-mz-to-db"]]);
const button = document.querySelector('.btn-grid');
const container = document.getElementById('add-to-database-adducts');

let divCounter = 0;

button?.addEventListener('click', () => {
  divCounter++;
  // Create new div
  const newDiv = document.createElement('div');
  newDiv.className = 'grid grid-4';

  // Add input fields to the div
  const input1 = document.createElement('input');
  input1.type = 'text';
  input1.placeholder = 'Adduct Name';
  input1.required = true
  input1.id = `add-adduct-name-${divCounter}`
  newDiv.appendChild(input1);

  const input2 = document.createElement('input');
  input2.type = 'text';
  input2.placeholder = 'Formula';
  input2.required = true;
  input2.id = `add-matrix-formula-${divCounter}`
  newDiv.appendChild(input2);

  const input3 = document.createElement('input');
  input3.type = 'text';
  input3.placeholder = 'No of derivs';
  input3.required = true;
  input3.id = `add-fg-${divCounter}`
  newDiv.appendChild(input3);

  // Add the reset button
  const resetButton = document.createElement('button');
  resetButton.type = 'reset';
  const spanBack = document.createElement('span');
  spanBack.className = 'back';
  const img = document.createElement('img');
  img.src = 'https://s3-us-west-2.amazonaws.com/s.cdpn.io/162656/eraser-icon.svg';
  spanBack.appendChild(img);
  resetButton.appendChild(spanBack);
  newDiv.appendChild(resetButton);

  resetButton.addEventListener('click', function() {
    newDiv.remove();
  });



  // Append new div to container
  container?.appendChild(newDiv);
});

document.getElementById('getValues')?.addEventListener('click', function() {
  const inputs = document.querySelectorAll('#container .grid input');
  const values: string[] = [];

  inputs.forEach((input) => {
    values.push((<HTMLInputElement>input).value);
  });

  console.log(values);  // Outputs the values to the console. Adjust this to process the values as needed.
});

if (slideInDiv?.classList.contains('slide-in--active')) {
    SubmitMatrixElement = document.getElementById("add-matrix-to-db") as HTMLButtonElement;
    MatrixRefreshElement = document.getElementById("ms1-slidein-table-h1-refresh-matrix") as HTMLSpanElement;
    createAndAttach("matrix");
    update_user_matrices();
}
}



export async function addToHTML_Matrix(adducts:string, table:string) {
  let checkbox_labels: string[] = await invoke("get_functional_groups_tauri", {});
		let checkboxes: string = generateCheckboxes(checkbox_labels, "Derivatizes:");
		adducts += 
		`
		<li class="add-to-database-adducts" id="add-to-database-adducts">
			<div class="grid grid-4">
				<input type="text" placeholder="Adduct Name" id="add-adduct-name-0" required>
				<input type="text" placeholder="Formula" id="add-matrix-formula-0" required>
				<input type="text" placeholder="No of derivs" id="add-fg-0" required>
				<button type="reset" enabled>
					<span class="back">
						<img src="https://s3-us-west-2.amazonaws.com/s.cdpn.io/162656/eraser-icon.svg" alt="">
					</span>
				</button> 
			</div>
		</li>

		<li>
			<div class="grid grid-1">
				<button class="btn-grid" type="text" placeholder="add row" enabled>
					<span class="front">Add row</span>
				</button>
				</div>
		</li>
		`;

		table += `
			<div class="ms1-slidein-spreadsheet-div" id="ms1-slidein-spreadsheet-div-matrix">
				<h1> MS1 User Matrices 
					<span class="ms1-slidein-table-h1-refresh-icon" id="ms1-slidein-table-h1-refresh-matrix">
						<ion-icon name="refresh"></ion-icon>
					</span>
				</h1>
				<div class="ms1-slidein-spreadsheet-table">
					<table class="ms1-slidein-datatable" id="ms1-slidein-datatable">
						<thead>
							<tr>
								<th>ID</th>
								<th>Name</th>
								<th>Derivatizes</th>
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
  return [checkboxes, adducts, table];
}




interface USER_MATRICES {
	ids: string[],
	names: string[],
	derivatizes: string[]
}

export async function update_user_matrices() {
	let matrices: USER_MATRICES = await invoke("update_user_matrices_tauri");

	await renderSmallItemsFor_UserMatrices(matrices[0], matrices[1], matrices[2])
}


export async function remove_row_from_user_matrices(row_id: string) {
	const parts = row_id.split("-");
	const lastpart = parts[parts.length - 1];
	const rowid = parseInt(lastpart);
	const done: number = await invoke("remove_row_from_user_matrices_tauri", {rowid:rowid});
	return done;
}

async function renderSmallItemsFor_UserMatrices(ids: string[], names: string[], derivatizes: string[]) {
	console.log(ids, names, derivatizes);
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
																<td>${derivatizes[i]}</td>
																<td>
																		<span class="ms1-slidein-table-trashbutton" id="ms1-slidein-table-trashbutton-${ids[i]}">
																				<ion-icon name="trash-outline"></ion-icon>
																		</span>
																</td>
															</tr>`;
				const templateContent = template.content;
				res_div.appendChild(templateContent);
                rowIds.push(`ms1-slidein-table-trashbutton-${ids[i]}`)
			}
        i++
		}

        removeFromDBInstance.addListeners(rowIds, "matrices")
		
		if (i < names.length-1) {
			setTimeout(appendBatch, 0);
		}
	}
	appendBatch();
}