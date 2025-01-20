import { listen } from "@tauri-apps/api/event";
import { addToHTML, add_to_db_rust, add_warning_message, cancellationToken, change_slidein, createAndAttach, generateCheckboxes, removeFromDBInstance, remove_warning_message, slideInDiv } from "../ms1_add_buttons";
import { invoke } from "@tauri-apps/api/core";


let FunctionalGroupListener: ((e: MouseEvent) => void) | null = null;
let FunctionalGroupRefreshListener: ((e: MouseEvent) => void) | null = null;



let SubmitFgElement: HTMLButtonElement | null = null;
let FgRefreshElement: HTMLSpanElement | null = null;


listen('update-progress', (event: any) => {
	const progress: number = event.payload;
	const progressBarFill: HTMLElement | null = document.getElementById('round-progress-bar-add-fg');
	if (progressBarFill) {
			const progressString = progress.toString();//.toFixed(2);
			progressBarFill.setAttribute('data-percent', progressString);
			updateProgressBarFromDataPercent(); // Refresh the progress bar visually

			if (progressString == "100")  {
				let addbutton = document.getElementById("add-fg-to-db") as HTMLButtonElement;
				addbutton.disabled = false;
				console.log("Button re-enabled");
			}
	}
});


function updateProgressBarFromDataPercent(): void {
	const el: HTMLElement | null = document.getElementById('round-progress-bar-add-fg');
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


export function createAndAttachFg() {
	FunctionalGroupListener = async () => {
		let name = document.getElementById("add-name-to-db") as HTMLInputElement;
		let smarts = document.getElementById("add-smarts-to-db") as HTMLInputElement;
		console.log("name:", name!.value);

		let duplicate = await invoke("check_fg_duplicate_tauri", {name:name!.value});
		console.log("duplicate: ", duplicate);
		if (duplicate == true) {
			console.log("found a duplicate")
			alert(`${name!.value} already exists`)
			return;
		} 
	
		const smartsResult = true; // Assume validation is successful for simplicity
	
		const fgCheckboxObject: Record<string, boolean> = {};
		const divElement = document.getElementById('add-to-database-checkboxes');
	
		if (divElement) {
			const checkboxes = divElement.querySelectorAll('input[type=checkbox]');
			checkboxes.forEach((checkbox) => {
				const inputCheckbox = checkbox as HTMLInputElement;
				fgCheckboxObject[inputCheckbox.id.replace("add-", "")] = inputCheckbox.checked;
			});
		}
	
		if (smartsResult) {
			remove_warning_message();
	
			try {
				// Disable the submit button
				let addbutton = document.getElementById("add-fg-to-db") as HTMLButtonElement;
				addbutton.disabled = true;
				console.log("Button disabled");
	
				// Wait for the Rust function to complete
				console.log("Starting Rust computation...");
				const result = await add_to_db_rust(
					name!.value,
					smarts!.value,
					"fg",
					fgCheckboxObject,
					{},
					[]
				);
				console.log("Rust computation complete:", result);
			} catch (error) {
				console.error("Error during computation:", error);
			} finally {
				// Re-enable the submit button
				//let addbutton = document.getElementById("add-fg-to-db") as HTMLButtonElement;
				//addbutton.disabled = false;
				//console.log("Button re-enabled");
			}
		} else {
			add_warning_message("SMARTS");
		}
	};
	
	SubmitFgElement?.addEventListener("click", FunctionalGroupListener);

	FunctionalGroupRefreshListener = async () => {
	update_user_fgs()
	}
	FgRefreshElement?.addEventListener("click", FunctionalGroupRefreshListener);
}

export function destroyAndDetachFg() {
	if (FunctionalGroupListener) {
		SubmitFgElement!.removeEventListener("click", FunctionalGroupListener);
		FunctionalGroupListener = null;
	}
	if (FunctionalGroupRefreshListener) {
		FgRefreshElement!.removeEventListener("click", FunctionalGroupRefreshListener);
		FunctionalGroupRefreshListener = null;
	}
}



export async function add_functional_group() {
	change_slidein();
	slideInDiv!.innerHTML = "";
	slideInDiv!.innerHTML = await addToHTML("add-fg-to-db", "Functional Group", [["Name", "add-name-to-db"], ["SMARTS", "add-smarts-to-db"]]);
	make_progressbar();

	if (slideInDiv?.classList.contains('slide-in--active')) {
		SubmitFgElement = document.getElementById("add-fg-to-db") as HTMLButtonElement;
		FgRefreshElement = document.getElementById("ms1-slidein-table-h1-refresh-fg") as HTMLSpanElement;
		createAndAttach("fg");
		update_user_fgs();
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



export async function addToHTML_FG(progressbar:string, table:string) {
  let checkbox_labels: string[] = [];
		checkbox_labels = await invoke("get_matrices_tauri", {});
		let checkboxes:string = generateCheckboxes(checkbox_labels, "Derivatized by:");
		progressbar += `
		<li>
			<div class="round-progress-bar" id="round-progress-bar-add-fg" data-percent="0"></div>
		</li>`;
		table += `
			<div class="ms1-slidein-spreadsheet-div" id="ms1-slidein-spreadsheet-div-fg">
				<h1> User functional groups 
					<span class="ms1-slidein-table-h1-refresh-icon" id="ms1-slidein-table-h1-refresh-fg">
						<ion-icon name="refresh"></ion-icon>
					</span>
				</h1>
				<div class="ms1-slidein-spreadsheet-table">
					<table class="ms1-slidein-datatable" id="ms1-slidein-datatable">
						<thead>
							<tr>
								<th>ID</th>
								<th>Name</th>
								<th>SMARTS</th>
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
      return [checkboxes, progressbar, table]
}





interface USER_FGS {
	ids: string[],
	names: string[],
	smarts: string[],
	finished: string[]
}

export async function update_user_fgs() {
	let fgs: USER_FGS = await invoke("update_user_fgs_tauri");

	await renderSmallItemsFor_UserFgs(fgs[0], fgs[1], fgs[2])
}


export async function remove_row_from_user_fgs(row_id: string, toremove:string) {
	const parts = row_id.split("-");
	const lastpart = parts[parts.length - 1];
	const rowid = parseInt(lastpart);
	const done: number = await invoke("remove_from_user_fgs_tauri", {rowid:rowid, toremove:toremove});
	return done;
}


async function renderSmallItemsFor_UserFgs(ids: string[], names: string[], smarts: string[]) {
	console.log("rendering small fgs");
	console.log(ids, names, smarts);
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
										<td>${smarts[i]}</td>
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

        removeFromDBInstance.addListeners(rowIds, "fgs")
		
		if (i < names.length-1) {
			setTimeout(appendBatch, 0);
		}
	}
	appendBatch();
}