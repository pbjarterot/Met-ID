import { invoke } from '@tauri-apps/api';
import { listen } from '@tauri-apps/api/event';

listen('update-progress', (event: any) => {
	const progress: number = event.payload;
	//console.log(progress);

	const progressBarFill: HTMLElement | null = document.getElementById('round-progress-bar-add-fg');
	if (progressBarFill) {
			const progressString = progress.toString();//.toFixed(2);
			console.log(progressString);
			progressBarFill.setAttribute('data-percent', progressString);
			updateProgressBarFromDataPercent(); // Refresh the progress bar visually
	}
});


function updateProgressBarFromDataPercent(): void {
	const el: HTMLElement | null = document.getElementById('round-progress-bar-add-fg');
	console.log(el);
	if (el) {
			renderProgressBar(el);
	}
}

function renderProgressBar(el: HTMLElement): void {
	const options = {
			percent: parseFloat(el.getAttribute('data-percent') || "25"), // Get it as a floating point number
			size: parseInt(el.getAttribute('data-size') || "220", 10),
			lineWidth: parseInt(el.getAttribute('data-line') || "15", 10),
			rotate: parseInt(el.getAttribute('data-rotate') || "0", 10)
	};

	console.log(options.percent);

	const canvas: HTMLCanvasElement = document.createElement('canvas');
	canvas.id = 'yourCanvasId';
	canvas.className = 'yourCanvasClass';

	const span: HTMLElement = document.createElement('span');
	span.id = 'yourSpanId';
	span.className = 'yourSpanClass';
	span.textContent = `${options.percent.toFixed(2)}%`;  // Use two decimal places for display

	const ctx = canvas.getContext('2d');
	if (ctx) {
			canvas.width = canvas.height = options.size;

			el.innerHTML = ""; // Clear previous content
			el.appendChild(span);
			el.appendChild(canvas);

			ctx.translate(options.size / 2, options.size / 2); // change center
			ctx.rotate((-1 / 2 + options.rotate / 180) * Math.PI); // rotate -90 deg

			const radius = (options.size - options.lineWidth) / 2;

			const drawCircle = (color: string, lineWidth: number, percent: number): void => {
					percent = Math.min(Math.max(0, percent), 1);
					ctx.beginPath();
					ctx.arc(0, 0, radius, 0, Math.PI * 2 * percent, false);
					ctx.strokeStyle = color;
					ctx.lineCap = 'round';
					ctx.lineWidth = lineWidth;
					ctx.stroke();
			};

			drawCircle('#5C5C5C', options.lineWidth, 100 / 100);
			drawCircle('#292929', options.lineWidth, options.percent / 100); // Round to whole number for drawing
	}
}



let MatrixListener: ((e: MouseEvent) => void) | null = null;
let MetaboliteListener: ((e: MouseEvent) => void) | null = null;
let FunctionalGroupListener: ((e: MouseEvent) => void) | null = null;

let SubmitMatrixElement: HTMLButtonElement | null = null;
let SubmitMetaboliteElement: HTMLButtonElement | null = null;
let SubmitFgElement: HTMLButtonElement | null = null;


const slideInDiv = document.getElementById('add-to-db-slidein');

function isNumeric(str: string): boolean {
	return !isNaN(Number(str));
}



async function add_to_db_rust(name: string, smilesSmartsMz: string, metType: string, endoExoOrOther, inTissue, adducts) {

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

function change_slidein() {
    if (slideInDiv!.classList.value === "slide-in") {
        slideInDiv!.classList.add('slide-in--active');
    } else if (slideInDiv!.classList.value === "slide-in slide-in--active"){
        slideInDiv!.classList.remove('slide-in--active');
		destroyAndDetach();
    }
}

function add_warning_message(type: string) {
	const parsingErrorDiv = document.getElementById("parsing-error") as HTMLDivElement;
	parsingErrorDiv!.textContent = type + " could not be parsed"; 
	parsingErrorDiv!.style.color = "red";
}

function remove_warning_message() {
	const parsingErrorDiv = document.getElementById("parsing-error") as HTMLDivElement;
	parsingErrorDiv!.textContent = "";

}

const createAndAttach = async (which_button: string) => {
	if (which_button === "matrix") {
		MatrixListener = async () => {
			let name = document.getElementById("add-name-to-db") as HTMLInputElement;
			let smiles_mz = document.getElementById("add-smiles-mz-to-db") as HTMLInputElement;
			const smilesResult = await check_smiles(smiles_mz!.value);

			const matrixCheckboxObject = {}

			const divElement = document.getElementById('add-to-database-checkboxes');
			if (divElement) {
				const checkboxes = divElement.querySelectorAll('input[type=checkbox]');
				checkboxes.forEach((checkbox) => {
					const inputCheckbox = checkbox as HTMLInputElement;
					matrixCheckboxObject[inputCheckbox.id.replace("add-", "")] = inputCheckbox.checked
					console.log(inputCheckbox.id, inputCheckbox.checked); // will be true or false
				});
			}

			const inputs = document.querySelectorAll('.add-to-database-adducts input');
			console.log(inputs);
			const values: string[] = [];

			inputs.forEach((input) => {
				values.push((<HTMLInputElement>input).value);
			});

			console.log(values);

			const matrixAdductsObject = {}

			for (let i = 0; i < inputs.length; i++) {
				let i_obj = inputs[i].id
				matrixAdductsObject[i_obj] = values[i];
			}
			console.log(matrixAdductsObject);


			if (smilesResult) {
				remove_warning_message();
				add_to_db_rust(name!.value, smiles_mz!.value, "matrix", matrixCheckboxObject, {}, matrixAdductsObject)
			} else if (isNumeric(smiles_mz!.value)) {
				remove_warning_message();
				add_to_db_rust(name!.value, smiles_mz!.value, "matrix", matrixCheckboxObject, {}, matrixAdductsObject)
			} else {
				add_warning_message("Value is not a number and SMILES")
			}


		};
		SubmitMatrixElement!.addEventListener("click", MatrixListener);
	}
	
	if (which_button === "metabolite") {
		MetaboliteListener = async () => {
			let name = document.getElementById("add-name-to-db") as HTMLInputElement;
			let smiles = document.getElementById("add-smiles-to-db") as HTMLInputElement;
			//const smilesResult = await check_smiles(smiles!.value);
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
					console.log(inputCheckbox.id, inputCheckbox.checked); // will be true or false
				});
			}
		
			if (smilesResult) {
				remove_warning_message();
				add_to_db_rust(name!.value, smiles!.value, "metabolite", MetTypeCheckboxObject, TissueCheckboxObject, {})
			} else {
				add_warning_message("SMILES");
			}

		};
		SubmitMetaboliteElement?.addEventListener("click", MetaboliteListener);
	}
	
	if (which_button === "fg") {
		FunctionalGroupListener = async () => {
			let name = document.getElementById("add-name-to-db") as HTMLInputElement;
			let smarts = document.getElementById("add-smarts-to-db") as HTMLInputElement;
			//const smartsResult = await check_smarts(smarts!.value);
			const smartsResult = true;

			const fgCheckboxObject = {}

			const divElement = document.getElementById('add-to-database-checkboxes');
			if (divElement) {
				const checkboxes = divElement.querySelectorAll('input[type=checkbox]');
				checkboxes.forEach((checkbox) => {
					const inputCheckbox = checkbox as HTMLInputElement;
					fgCheckboxObject[inputCheckbox.id.replace("add-", "")] = inputCheckbox.checked
					console.log(inputCheckbox.id, inputCheckbox.checked); // will be true or false
				});
			}

			if (smartsResult) {
				remove_warning_message();
				add_to_db_rust(name!.value, smarts!.value, "fg", fgCheckboxObject, {}, {})
			} else {
				add_warning_message("SMARTS")
			}

		};
		SubmitFgElement?.addEventListener("click", FunctionalGroupListener);
	}
};

const destroyAndDetach = () => {
	// Remove event listeners
	if (MatrixListener) {
		SubmitMatrixElement!.removeEventListener("click", MatrixListener);
		MatrixListener = null;
	}
	if (MetaboliteListener) {
		SubmitMetaboliteElement!.removeEventListener("click", MetaboliteListener);
		MetaboliteListener = null;
	}
	if (FunctionalGroupListener) {
		SubmitFgElement!.removeEventListener("click", FunctionalGroupListener);
		FunctionalGroupListener = null;
	}
}

async function check_smiles(smiles: string): Promise<boolean> {
    let smiles_ok: boolean = await invoke("check_smiles", {smiles:smiles});

    return smiles_ok;
}
/*
async function check_smarts(smarts: string): Promise<boolean> {
    let smarts_ok: boolean = await invoke("check_smarts", {smarts:smarts});

    return smarts_ok;
}
*/

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
    	createAndAttach("matrix");
	}
}

export async function add_metabolite() {
    change_slidein();
	slideInDiv!.innerHTML = "";
	slideInDiv!.innerHTML = await addToHTML("add-metabolite-to-db", "Metabolite", [["Name", "add-name-to-db"], ["SMILES", "add-smiles-to-db"]]);

	if (slideInDiv?.classList.contains('slide-in--active')) {
		SubmitMetaboliteElement = document.getElementById("add-metabolite-to-db") as HTMLButtonElement;
		createAndAttach("metabolite");
	}
}

export async function add_functional_group() {
	change_slidein();
	slideInDiv!.innerHTML = "";
	slideInDiv!.innerHTML = await addToHTML("add-fg-to-db", "Functional Group", [["Name", "add-name-to-db"], ["SMARTS", "add-smarts-to-db"]]);
	make_progressbar();

	if (slideInDiv?.classList.contains('slide-in--active')) {
		SubmitFgElement = document.getElementById("add-fg-to-db") as HTMLButtonElement;
		createAndAttach("fg");
	}
}

function make_progressbar() {
	const el: HTMLElement | null = document.getElementById('round-progress-bar-add-fg'); // get canvas
	console.log(el);

	if (el) {
			const options = {
					percent: parseInt(el.getAttribute('data-percent') || "25", 10),
					size: parseInt(el.getAttribute('data-size') || "220", 10),
					lineWidth: parseInt(el.getAttribute('data-line') || "15", 10),
					rotate: parseInt(el.getAttribute('data-rotate') || "0", 10)
			};

			const canvas: HTMLCanvasElement = document.createElement('canvas');
			canvas.id = 'yourCanvasId';
			canvas.className = 'yourCanvasClass';
			const span: HTMLElement = document.createElement('span');
			span.id = 'yourSpanId';
			span.className = 'yourSpanClass';
			span.textContent = `${options.percent}%`;

			const ctx = canvas.getContext('2d');
			if (ctx) {
					canvas.width = canvas.height = options.size;

					el.appendChild(span);
					el.appendChild(canvas);

					ctx.translate(options.size / 2, options.size / 2); // change center
					ctx.rotate((-1 / 2 + options.rotate / 180) * Math.PI); // rotate -90 deg

					const radius = (options.size - options.lineWidth) / 2;

					const drawCircle = (color: string, lineWidth: number, percent: number): void => {
							percent = Math.min(Math.max(0, percent), 1);
							ctx.beginPath();
							ctx.arc(0, 0, radius, 0, Math.PI * 2 * percent, false);
							ctx.strokeStyle = color;
							ctx.lineCap = 'round'; // butt, round or square
							ctx.lineWidth = lineWidth;
							ctx.stroke();
					};

					drawCircle('#5C5C5C', options.lineWidth, 100 / 100);
					drawCircle('#292929', options.lineWidth, options.percent / 100);
			}
	}
}


async function addToHTML(submit_id: string, add_term: string, inputs: string[][]) {
	
	let HTMLInputs = "";
	let progressbar = "";
	//let adducts = "";
	for (let i = 0; i < inputs.length; i++) {
		HTMLInputs += `<input type="text" placeholder="${inputs[i][0]}" id="${inputs[i][1]}" required>`
	}

	let checkboxes = "";
	let adducts = "";
	if (add_term === "Metabolite") {
		checkboxes = generateCheckboxes(["Endogenous", "Exogenous", "Unspecified"], "Type of Metabolite");

		let checkbox_labels: string[] = await invoke("get_tissues", {});
		checkboxes += generateCheckboxes(checkbox_labels, "Found in tissue:");


	} else if (add_term === "Matrix") {
		let checkbox_labels: string[] = await invoke("get_functional_groups", {});
		checkboxes = generateCheckboxes(checkbox_labels, "Derivatizes:");
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

	} else if (add_term === "Functional Group") {
		let checkbox_labels: string[] = [];
		checkbox_labels = await invoke("get_matrices", {});
		checkboxes = generateCheckboxes(checkbox_labels, "Derivatized by:");
		progressbar += `
		<li>
			<div class="round-progress-bar" id="round-progress-bar-add-fg" data-percent="0"></div>
		</li>`;
	} else {
		//return ""
		console.log("something bad happened");
	}
	
	
	let innerHTML = `
	<div class="add-to-database">
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
	
	`;


	
	return innerHTML;
}

function generateCheckboxes(labels: string[], title: string): string {
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