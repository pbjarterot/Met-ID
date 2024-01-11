import { invoke } from "@tauri-apps/api/tauri";
import { open } from '@tauri-apps/api/dialog';
import { generate_ms2_small_results_card, small_results_card_expanding } from "./ms2_small_results_card";
import { compare_msms } from "./ms2_compare";
import { get_msms_from_mzml, match_msms } from "./ms2_io";
import { showPopup } from "./ms2_popup";

//var stopRendering: boolean = false;

var slider2 = document.getElementById("mzWindow2") as HTMLInputElement;
var output2 = document.getElementById("demo2") as HTMLElement;
    output2!.innerHTML = slider2!.value + " mDa";

slider2!.oninput = function() {
    output2!.innerHTML = (this as HTMLInputElement).value + " mDa";
}

var slider3 = document.getElementById("mzWindow3") as HTMLInputElement;
var output3 = document.getElementById("demo3") as HTMLElement;
    output3!.innerHTML = slider3!.value + " mDa";

slider3!.oninput = function() {
    output3!.innerHTML = (this as HTMLInputElement).value + " mDa";
}

export interface MSMSDatabase {
	names: string[],
	identifiers: string[],
	adducts: string[],
	cids: string[],
	windows: string[],
	tofs: string[],
	mzs: string[],
	cossim: number[],
	matrices: string[]
}


const checkBackendReadiness = async () => {
  try {
    const isReady = await invoke('is_backend_ready');
    if (isReady) {
      // Backend is ready, proceed with your app logic.
    } else {
      // Wait for a short period and check again.
      setTimeout(checkBackendReadiness, 10000);
    }
  } catch (error) {
    console.log('Error checking backend readiness:', error);
  }
}

function delay(ms: number): Promise<void> {
	return new Promise(resolve => setTimeout(resolve, ms));
}

window.addEventListener("DOMContentLoaded", async () => {
	checkBackendReadiness();
	
	document.getElementById("ms2-sidebar-open-file-button")!.addEventListener("click", async () => get_msms_from_mzml());
	document.getElementById("ms2-sidebar-match-button")!.addEventListener("click", async () => match_msms());
	document.getElementById("ms2-sidebar-add-to-db-button")!.addEventListener("click", async () => add_msms_to_db());
	
	await delay(100)
	let [names, identifiers, adducts, cids, mzs, _windows, _tofs, matrices]: [string[], string[], string[], string[], string[], string[], string[], string[]] = await invoke("get_msms_tauri", {});

	let msms: MSMSDatabase ={
		names: names,
		identifiers: identifiers,
		adducts: adducts,
		cids: cids,
		windows: _windows,
		tofs: _tofs,
		mzs: mzs,
		cossim: new Array(mzs.length).fill(0),
		matrices: matrices
	};

	await renderSmallItems(msms.names, msms.cids, msms.adducts, msms.mzs, msms.identifiers, msms.cossim, msms.matrices);
	document.getElementById("ms2-compare")!.addEventListener("click", () => compare_msms());

	//Searchbar listeners
	addMSMSSearchbarListeners();

});


async function add_msms_to_db() {
	let result = await open({ directory: false, multiple: false, filters: [{ name: 'mzML', extensions: ['mzML'] }] }) as unknown as string;

	showPopup(result);

}

function addMSMSSearchbarListeners() {
	// Fetch the parent div by its id
	const parentDiv = document.getElementById("ms2-searchbar-left");

	if (parentDiv) {
		// Listen for input events on each child element within the div
		parentDiv.addEventListener("input", (event) => {
		// Optional: Check if the event's target is an input element
		if (event.target instanceof HTMLInputElement) {
			// Trigger your update function
			//stopRendering = true;
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
		});
	}
}

function isEmpty(obj: object): boolean {
    return Object.keys(obj).length === 0;
}

// Your update function
export async function updateMSMSResults(sorted: boolean, msms: MSMSDatabase) {
	cancellationToken.isCancelled = false;

	
	let name = (document.getElementById("msms-searchbar") as HTMLInputElement)!.value;
	let fragment = (document.getElementById("msms-searchbar-fragment") as HTMLInputElement)!.value;
	let ms1mass = (document.getElementById("msms-searchbar-ms1") as HTMLInputElement)!.value;

	let fragmentslider = (document.getElementById("mzWindow2") as HTMLInputElement)!.value;
	let ms1massslider = (document.getElementById("mzWindow3") as HTMLInputElement)!.value;

	if (cancellationToken.isCancelled) {
		console.log("Operation cancelled");
		return;
	}
	if (sorted === false) {
		let [names, identifiers, adducts, cids, windows, tofs, mzs, _cossim, matrices]: [string[], string[], string[], string[], string[], string[], string[], number[], string[]] = await invoke("ms2_search_spectra_tauri", {name: name, fragment:fragment, ms1mass:ms1mass, fragmentslider:fragmentslider, ms1massslider:ms1massslider});
		let a: MSMSDatabase ={
			names: names,
			identifiers: identifiers,
			adducts: adducts,
			cids: cids,
			windows: windows,
			tofs: tofs,
			mzs: mzs,
			cossim: new Array(mzs.length).fill(0),
			matrices: matrices
		};
		msms = a;
	} 

	// Check the cancellation token again
	if (cancellationToken.isCancelled) {
		console.log("Operation cancelled");
		return;
	}

	if (isEmpty(msms)) {
		console.log("msms empty")
		return;
	}
													//names			cids			adducts					mzs			identifiers				cossim
	await renderSmallItems(msms.names, msms.cids, msms.adducts, msms.mzs, msms.identifiers, msms.cossim, msms.matrices);
}

// Create a cancellation token object
export interface CancellationToken {
	isCancelled: boolean;
}

let cancellationToken: CancellationToken = { isCancelled: false };


async function renderSmallItems(names: string[], cids: string[], adducts: string[], parent_mzs: string[] | number[], identifiers: string[], cossim: number[], matrices: string[]) {
	const res_div = document.querySelector("#ms2-results-body") as HTMLDivElement;
	res_div.innerHTML = "";
	const batchSize = 3;
	let i = 0;

	type myTuple = [string, string, string];
	let already_done: myTuple[] = [];


	function appendBatch() {
		if (cancellationToken.isCancelled) {
			console.log("Operation cancelled");
			return;
		}
		const template = document.createElement("template");
		for (let j = 0; j < batchSize; j++) {
			//checking for 0eV is not optimal when searching, as fragments could only appear in e.g. 30eV
			if (cancellationToken.isCancelled) {
				console.log("Operation cancelled");
				return;

			}
			if (!already_done.some(([name, adduct, cid]) => name === names[i] && adduct === adducts[i] && cid === cids[i])) {
				template.innerHTML = generate_ms2_small_results_card(names[i], cids[i], adducts[i], parent_mzs[i], cossim[i], matrices[i], i);
				const templateContent = template.content;
				res_div.appendChild(templateContent);
				let a = document.getElementById("ms2-small-results-card-" + i.toString()) as HTMLDivElement;
				let b = a.querySelector("div.ms2-small-results-card-name") as HTMLDivElement;
				small_results_card_expanding(a, b, identifiers[i], i, adducts[i], cids[i]);
				already_done.push([names[i], adducts[i], cids[i]])
			}
        i++
		}
		
		if (i < names.length-1) {
			setTimeout(appendBatch, 0);
		}
	}
	appendBatch();
}


//document.getElementById("ms2-sidebar-open-file-button")?.addEventListener("click", () => openFile());


export interface All_spectra_database {
	name: string[],
	//matrix: string[],
	smiles: string[],
	cid: string[],
	fragment: string[]
}

