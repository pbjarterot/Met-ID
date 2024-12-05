import { invoke } from '@tauri-apps/api/core';
import { CancellationToken } from '../ms2/ms2_main';
import { addToHTML_Metabolite, createAndAttachMetabolites, destroyAndDetachMetabolite, remove_row_from_user_metabolites, update_user_metabolites } from './add_buttons/metabolites';
import { addToHTML_Matrix, createAndAttachMatrices, destroyAndDetachMatrix, remove_row_from_user_matrices, update_user_matrices } from './add_buttons/matrices';
import { addToHTML_FG, createAndAttachFg, destroyAndDetachFg, remove_row_from_user_fgs, update_user_fgs } from './add_buttons/functional_groups';
import { fill_dropdown, new_tgt_matrix } from '../dropdown';

export let cancellationToken: CancellationToken = { isCancelled: false };

export const slideInDiv = document.getElementById('add-to-db-slidein');



export async function add_to_db_rust(name: string, smilesSmartsMz: string, metType: string, endoExoOrOther:object, inTissue:object, adducts:string[]) {

	let added_to_db = await invoke("my_command", {
		args: { 
			name: name, 
			smilesSmartsMz: smilesSmartsMz, 
			metType:metType, 
			endoExoOrOther: endoExoOrOther, 
			inTissue: inTissue,
			adducts: adducts
		}
	});
	return added_to_db;
}

export function change_slidein() {
    if (slideInDiv!.classList.value === "slide-in") {
        slideInDiv!.classList.add('slide-in--active');
    } else if (slideInDiv!.classList.value === "slide-in slide-in--active"){
        slideInDiv!.classList.remove('slide-in--active');
		destroyAndDetach();
    }
}

export function add_warning_message(type: string) {
	const parsingErrorDiv = document.getElementById("parsing-error") as HTMLDivElement;
	parsingErrorDiv!.textContent = type; 
	parsingErrorDiv!.style.color = "red";
}

export function remove_warning_message() {
	const parsingErrorDiv = document.getElementById("parsing-error") as HTMLDivElement;
	parsingErrorDiv!.textContent = "";

}

export const createAndAttach = async (which_button: string) => {
	if (which_button === "matrix") {
		createAndAttachMatrices();
	}
	if (which_button === "metabolite") {
		createAndAttachMetabolites();
	}
	if (which_button === "fg") {
		createAndAttachFg();
	}
};

const destroyAndDetach = () => {
	// Remove event listeners
	destroyAndDetachMatrix();

	destroyAndDetachMetabolite();

	destroyAndDetachFg();
	
}



export async function addToHTML(submit_id: string, add_term: string, inputs: string[][]) {
	
	let HTMLInputs = "";
	let progressbar = "";
	let table = "";
	//let adducts = "";
	for (let i = 0; i < inputs.length; i++) {
		HTMLInputs += `<input type="text" placeholder="${inputs[i][0]}" id="${inputs[i][1]}" required>`
	}

	let checkboxes = "";
	let adducts = "";
	if (add_term === "Metabolite") {
		[checkboxes, table] = await addToHTML_Metabolite(table);

	} else if (add_term === "Matrix") {
		[checkboxes, adducts, table] = await addToHTML_Matrix(adducts, table);

	} else if (add_term === "Functional Group") {
		[checkboxes, progressbar, table] = await addToHTML_FG(progressbar, table);
	} else {
		//return ""
	}
	
	let innerHTML = `
	<div class="add-to-database">
		<div class="add-to-database-top">
			<h1>Add ${add_term}</h1>
			<li>
				<div class="grid grid-1">
					${HTMLInputs}
				</div>
			</li>
			<li>
				<!--
					<div class="required-msg">ALL FIELDS ARE REQUIRED</div>
				-->
				<div class="parsing-err" id="parsing-error"></div>

			</li>
			<li class="add-to-database-checkboxes" id="add-to-database-checkboxes">
				${checkboxes}
				
			</li>

			${adducts}
			
			<li>
				<div class="grid grid-2">
					<button class="btn-grid" type="submit" id="${submit_id}" enabled>
						<span class="back">
							<img src="https://s3-us-west-2.amazonaws.com/s.cdpn.io/162656/email-icon.svg" alt="">
						</span>
						<span class="front">SUBMIT</span>
					</button>
					<button class="btn-grid" type="reset" disabled>
						<span class="back">
							<img src="https://s3-us-west-2.amazonaws.com/s.cdpn.io/162656/eraser-icon.svg" alt="">
						</span>
						<span class="front">RESET</span>
					</button> 
				</div>
			</li>  
			
			${progressbar}
		</div>
		${table}
	</div>
	
	`;

	return innerHTML;
}

export function generateCheckboxes(labels: string[], title: string): string {
	let checkboxesHtml = `<h2>${title}</h2><div class="grid">`;

	// Function to measure the width of text
	const measureText = (text: string, fontSize: number) => {
		const canvas = document.createElement("canvas");
		const context = canvas.getContext("2d");
		if (context) {
			context.font = `${fontSize}px Arial`;
		return context.measureText(text).width;
		}
		return 0;
	};

	// You should replace this with the actual width of your container
	const containerWidth = 360;

	let remainingWidth = containerWidth;

	for (const label of labels) {
		// Assuming 12px font size and additional 20px for checkbox and gap
		const labelWidth = measureText(label, 12) + 30; 

		if (labelWidth > remainingWidth) {
		checkboxesHtml += `</div><div class="grid">`;
		remainingWidth = containerWidth;
		}

		checkboxesHtml += `
		<div class="add-to-database-checkbox-child">
		<input type="checkbox" id="add-${label}">
		<label for="add-${label}">${label}</label>
		</div>`;

		remainingWidth -= (labelWidth + 30);  // Assume 2px gap between checkboxes
	}

	checkboxesHtml += '</div>';
	return checkboxesHtml;
	}


class removeFromDBListeners {
	public removeListenerFunctions: Array<() => void> = [];
	public addListeners(elementIds: string[], mode: string) {
			const events = ['click']; // Example events
			console.log("elemids", elementIds)
			elementIds.forEach(id => {
			events.forEach(eventType => {
					const element = document.getElementById(id);
					if (!element) return;

					const listener = async () => {
							console.log(id)
							if (mode === "metabolites") {
								remove_row_from_user_metabolites(id);
								this.removeListeners();
								update_user_metabolites();
							} else if (mode === "matrices") {
								remove_row_from_user_matrices(id);
								this.removeListeners();
								update_user_matrices();
								fill_dropdown(Object.keys(await new_tgt_matrix()), "matrix-dropdown")
							} else if (mode === "fgs") {
								let name = get_name_in_table(element);
								if (name) {
									remove_row_from_user_fgs(id, name);
								}
								
								this.removeListeners();
								update_user_fgs();
							}
							
					}
					element.addEventListener(eventType, listener);
					this.removeListenerFunctions.push(() => element.removeEventListener(eventType, listener));
					});
			});
	}

	public removeListeners() {
			this.removeListenerFunctions.forEach(removeListener => removeListener());
			this.removeListenerFunctions = []; // Clear the array
	}
}

export const removeFromDBInstance = new removeFromDBListeners();



function get_name_in_table(childElement: HTMLElement) {

	let secondTDText: string | null = "";
	if (childElement) {
    // Getting the parent <tr> element
    let parentTR = childElement.closest('tr');

    if (parentTR) {
        // Getting all <td> elements in the parent <tr>
        let tds = parentTR.querySelectorAll('td');

        // Checking if there are at least two <td> elements
        if (tds.length >= 2) {
            // Getting the text content of the second <td>
            secondTDText = tds[1].textContent;
        }
    }
}
	return secondTDText;
}