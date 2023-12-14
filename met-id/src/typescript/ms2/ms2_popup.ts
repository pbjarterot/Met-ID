let AddToDbListener: ((e: MouseEvent) => void) | null = null;
import { invoke } from '@tauri-apps/api/tauri';

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
        }
    }

    const addToDbButton = document.getElementById("ms2-add-to-db-add-button");
    if (AddToDbListener) {
        addToDbButton!.addEventListener("click", AddToDbListener);
      }

}

function hidePopup() {
	const overlay = document.getElementById('ms2-add-to-db-overlay');
	const popup = document.getElementById('ms2-add-to-db-popup');

	overlay!.style.display = 'none';
	popup!.style.display = 'none';

    const addToDbButton = document.getElementById("ms2-add-to-db-add-button");

    if (AddToDbListener) {
        addToDbButton!.removeEventListener("click", AddToDbListener);
        AddToDbListener = null;
      }
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

    const _parsed_mzml: string = await invoke("add_msms_to_db", {name:name, adduct:adduct, mz:mz, cid:cid, tof:tof, mzwindow:mzwindow, identifier:identifier });
    
    console.log("adding", _parsed_mzml);
    return
}


export function add_inputs_to_popup() {
	let inputdiv = document.getElementById("ms2-add-to-db-div");
	
	let HTMLInputs = "";
	HTMLInputs += `<div class="grid grid-2"><input type="text" placeholder="Name" id="ms2-popup-add-name-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="Identifier" id="ms2-popup-add-id-to-db" required></div>`;
	HTMLInputs += `<div class="grid grid-5"><input type="text" placeholder="Adduct" id="ms2-popup-add-adduct-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="cid" id="ms2-popup-add-cid-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="window" id="ms2-popup-add-window-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="tof" id="ms2-popup-add-tof-to-db" required>`;
	HTMLInputs += `<input type="text" placeholder="mz" id="ms2-popup-add-mz-to-db" required></div>`;
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
    return 1;
}