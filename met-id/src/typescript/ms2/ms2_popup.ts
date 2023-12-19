
import { invoke } from '@tauri-apps/api/tauri';
import { CancellationToken, MSMSDatabase, updateMSMSResults } from './ms2_main';

let cancellationToken: CancellationToken = { isCancelled: false };
let AddToDbListener: ((e: MouseEvent) => void) | null = null;

let continuelistener: ((e: MouseEvent) => void) | null = null;
let cancellistener: ((e: MouseEvent) => void) | null = null;


class removeFromDBListeners {

    public removeListenerFunctions: Array<() => void> = [];


    public addListeners(elementIds: string[]) {
        const events = ['click']; // Example events
        //const elementIds = ['element1', 'element2', 'element3']; // Example element IDs

        elementIds.forEach(id => {
        events.forEach(eventType => {
            const element = document.getElementById(id);
            if (!element) return;

            const listener = () => {
                console.log(`Event ${eventType} on ${id}`);
                let done = remove_row_from_msms_db(id);
                console.log(done);
                this.removeListeners();
                update_user_msms();
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

const removeFromDBInstance = new removeFromDBListeners();

async function update_user_msms() {

    const msms: MSMSDatabase = await invoke("show_user_msms_db", {});
    
    console.log(msms);
    await renderSmallItems_for_user_msms(msms[0], msms[1], msms[2], msms[3], msms[4], msms[5], msms[6], msms[7], msms[8]);
}

async function remove_row_from_msms_db(row_id: string) {
    const parts = row_id.split("-");
    const lastpart = parts[parts.length - 1];
    const rowid = parseInt(lastpart);
    const done: number = await invoke("remove_row_from_msms_user_db", {rowid:rowid});
    return done;
}

async function renderSmallItems_for_user_msms(ids: string[], names: string[], identifiers: string[], adducts: string[], cids: string[] | number[], mzwindows: string[], tofs: string[], mzs: string[], matrices: string[]) {
	const res_div = document.getElementById("ms2-popup-tbody") as HTMLTableElement;
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
                                        <td>${identifiers[i]}</td>
                                        <td>${adducts[i]}</td>
                                        <td>${cids[i]}</td>
                                        <td>${mzwindows[i]}</td>
                                        <td>${tofs[i]}</td>
                                        <td>${mzs[i]}</td>
                                        <td>${matrices[i]}</td>
                                        <td>
                                            <button class="ms2-add-to-db-trashbutton" id="ms2-add-to-db-trashbutton-${ids[i]}">
                                                <ion-icon name="trash-outline"></ion-icon>
                                            </button>
                                        </td>
                                    </tr>`;
				const templateContent = template.content;
				res_div.appendChild(templateContent);
                rowIds.push(`ms2-add-to-db-trashbutton-${ids[i]}`)
			}
        i++
		}

        removeFromDBInstance.addListeners(rowIds)
		
		if (i < names.length-1) {
			setTimeout(appendBatch, 0);
		}
	}
	appendBatch();
}

export function showPopup (result: string) {
	console.log("showPopup");
    const overlay = document.getElementById('ms2-add-to-db-overlay');
    const popup = document.getElementById('ms2-add-to-db-popup');

    overlay!.style.display = 'block';
    popup!.style.display = 'flex';
    
    const closePopupButton = document.getElementById("ms1-popup-cancel-button");
    closePopupButton?.addEventListener("click", () => hidePopup())

    const overlayDiv = document.getElementById("ms2-add-to-db-overlay");
    overlayDiv?.addEventListener("click", () => hidePopup());

    add_filename_to_top_of_popup(result)
	add_inputs_to_popup();

    AddToDbListener = () => {
        if (checkInputs() === 0) {
            return;
        } else {
            add_to_msms_db();
            update_user_msms();
        }
    }

    continuelistener = () => {
        hidePopup();
        let a: MSMSDatabase = {
            names: [],
            identifiers: [],
            adducts: [],
            cids: [],
            windows: [],
            tofs: [],
            mzs: [],
            cossim: [],
            matrices: []
        };
        updateMSMSResults(false, a);
    }
    cancellistener = () => {
        hidePopup()
    }


    const continueButton = document.getElementById("ms2-popup-continue-button");
    if (continuelistener) {
        continueButton!.addEventListener("click", continuelistener);
    }
    const cancelButton = document.getElementById("ms2-popup-cancel-button")
    if (cancellistener) {
        cancelButton!.addEventListener("click", cancellistener);
    }

    const addToDbButton = document.getElementById("ms2-add-to-db-add-button");
    if (AddToDbListener) {
        addToDbButton!.addEventListener("click", AddToDbListener);
    }

    update_user_msms()

}

function hidePopup() {
	const overlay = document.getElementById('ms2-add-to-db-overlay');
	const popup = document.getElementById('ms2-add-to-db-popup');

	overlay!.style.display = 'none';
	popup!.style.display = 'none';

    const continueButton = document.getElementById("ms2-popup-continue-button");
    if (continuelistener) {
        continueButton!.removeEventListener("click", continuelistener);
    }
    const cancelButton = document.getElementById("ms2-popup-cancel-button")
    if (cancellistener) {
        cancelButton!.removeEventListener("click", cancellistener);
    }



    const addToDbButton = document.getElementById("ms2-add-to-db-add-button");

    if (AddToDbListener) {
        addToDbButton!.removeEventListener("click", AddToDbListener);
        AddToDbListener = null;
    }

    removeFromDBInstance.removeListeners()
}

function add_filename_to_top_of_popup(result: string){
	let resultString: string | null = Array.isArray(result) ? result.join(", ") : result;
	let filenameWindow = document.getElementById("ms2-popup-filename-window-name");
	filenameWindow!.textContent = resultString;
}

async function add_to_msms_db() {
    var name = (document.getElementById('ms2-popup-add-name-to-db')! as HTMLInputElement).value;
    var identifier = (document.getElementById('ms2-popup-add-id-to-db')! as HTMLInputElement).value;
    var adduct = (document.getElementById('ms2-popup-add-adduct-to-db')! as HTMLInputElement).value;
    var cid = (document.getElementById('ms2-popup-add-cid-to-db')! as HTMLInputElement).value;
    var mzwindow = (document.getElementById('ms2-popup-add-window-to-db')! as HTMLInputElement).value;
    var tof = (document.getElementById('ms2-popup-add-tof-to-db')! as HTMLInputElement).value;
    var mz = (document.getElementById('ms2-popup-add-mz-to-db')! as HTMLInputElement).value;
    var matrix = (document.getElementById('ms2-popup-add-matrix-to-db')! as HTMLInputElement).value;
    var path = document.getElementById('ms2-popup-filename-window-name')!.textContent;
    console.log(path);

    const _parsed_mzml: string = await invoke("add_msms_to_db", {name:name, adduct:adduct, mz:mz, cid:cid, tof:tof, mzwindow:mzwindow, identifier:identifier, path:path, matrix:matrix});
    
    console.log("adding", _parsed_mzml);
    return
}

export function add_inputs_to_popup() {
	let inputdiv = document.getElementById("ms2-add-to-db-div");
	
	let HTMLInputs = "";
	HTMLInputs += `<div class="grid grid-2"><input type="text" placeholder="Name" id="ms2-popup-add-name-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="Identifier" id="ms2-popup-add-id-to-db" required></div>`;
	HTMLInputs += `<div class="grid grid-6"><input type="text" placeholder="Adduct" id="ms2-popup-add-adduct-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="cid" id="ms2-popup-add-cid-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="window" id="ms2-popup-add-window-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="tof" id="ms2-popup-add-tof-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="mz" id="ms2-popup-add-mz-to-db" required>`;
    HTMLInputs += `<input type="text" placeholder="matrix" id="ms2-popup-add-matrix-to-db" required></div>`;
    HTMLInputs += `<button class="grid grid-1" id="ms2-add-to-db-add-button">Add to Database</button>`
	inputdiv!.innerHTML = HTMLInputs;
}

function checkInputs() {
    var input1 = (document.getElementById('ms2-popup-add-name-to-db')! as HTMLInputElement).value;
    var input2 = (document.getElementById('ms2-popup-add-id-to-db')! as HTMLInputElement).value;
    var input3 = (document.getElementById('ms2-popup-add-adduct-to-db')! as HTMLInputElement).value;
    var input4 = (document.getElementById('ms2-popup-add-cid-to-db')! as HTMLInputElement).value;
    var input5 = (document.getElementById('ms2-popup-add-window-to-db')! as HTMLInputElement).value;
    var input6 = (document.getElementById('ms2-popup-add-tof-to-db')! as HTMLInputElement).value;
    var input7 = (document.getElementById('ms2-popup-add-mz-to-db')! as HTMLInputElement).value;
    var input8 = (document.getElementById('ms2-popup-add-matrix-to-db')! as HTMLInputElement).value;
    
    if (!input1) {
        alert('Name is Empty');
        return 0;
    }
    if (!input2) {
        alert('ID is Empty');
        return 0;
    }
    if (!input3) {
        alert('Adduct is Empty');
        return 0;
    }
    if (!input4) {
        alert('CID is Empty');
        return 0;
    }
    if (!input5) {
        alert('Window is Empty');
        return 0;
    }
    if (!input6) {
        alert('TOF is Empty');
        return 0;
    }
    if (!input7) {
        alert('m/z is Empty');
        return 0;
    }
    if (!input8) {
        alert('Matrix is Empty');
        return 0;
    }
    return 1;
}