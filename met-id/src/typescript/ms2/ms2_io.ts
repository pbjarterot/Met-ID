import { open } from '@tauri-apps/api/dialog';
import { invoke } from '@tauri-apps/api/tauri';
//import { generate_ms2_small_results_card, small_results_card_expanding } from "./ms2_small_results_card";
import { MSMSDatabase, updateMSMSResults } from './ms2_main';

let result: string | string[] | null = "";


let MatchListener: ((e: MouseEvent) => void) | null = null;

async function read_input_mzml(result: string) {
  const parsed_mzml: string = await invoke("read_mzml_for_msms", {path: result});
  console.log("parsed_csv: ", parsed_mzml)
  return parsed_mzml
}

export async function get_msms_from_mzml() {
    result = await open({directory: false, multiple: false, filters: [{name: 'mzML', extensions:['mzML']}]}) as string;

    if (result && result.length > 0) {
      read_input_mzml(result);
      let a: MSMSDatabase ={
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

      //addUserSpectra()
    }
}
export function showMatchPopup (this: any) {
    const overlay = document.getElementById('ms2-match-overlay');
    const popup = document.getElementById('ms2-match-popup');

    overlay!.style.display = 'block';
    popup!.style.display = 'flex';
    
    //const closePopupButton = document.getElementById("ms1-popup-cancel-button");
    //closePopupButton?.addEventListener("click", () => hideMatchPopup())
    const popup_options = document.getElementById("ms2-match-options");

    popup_options!.innerHTML = `<div class="ms2-match-option-row">
                                  <div class="ms2-match-option-binsize">
                                    <p>Bin size (Da)</p>
                                  </div>
                                  <input type="text" class="ms2-match-option-input" id="ms2-match-option-binsize-input" oninput="validateInput(this)">
                                </div>`;
    const inputElement = document.getElementById("ms2-match-option-binsize-input");

    if (inputElement) {
        inputElement.oninput = function(event) {
            validateInput(event.target as HTMLInputElement);
        };
    } else {
        console.error("Element with ID 'ms2-match-option-binsize-input' not found.");
    }
                                

    const continuePopupButton = document.getElementById("ms2-match-popup-continue-button");
    if (MatchListener) {
      continuePopupButton!.addEventListener("click", MatchListener);
    }
    

    const overlayDiv = document.getElementById("ms2-match-overlay");
    overlayDiv?.addEventListener("click", () => hideMatchPopup());

}

function validateInput(input: HTMLInputElement): void {
  const value = input.value;
  const numbers = value.split('').filter(char => !isNaN(Number(char)) || char === '.').join('');
  let hasDecimal = false;
  let filtered = '';

  for (const char of numbers) {
      if (char === '.' && !hasDecimal) {
          hasDecimal = true;
          filtered += char;
      } else if (char !== '.') {
          filtered += char;
      }
  }

  input.value = filtered;
}


function hideMatchPopup() {
	const overlay = document.getElementById('ms2-match-overlay');
	const popup = document.getElementById('ms2-match-popup');

	overlay!.style.display = 'none';
	popup!.style.display = 'none';

  const continuePopupButton = document.getElementById("ms2-match-popup-continue-button");
  if (MatchListener) {
    continuePopupButton!.removeEventListener("click", MatchListener);
    MatchListener = null;
  }
    
}

async function run_msms_match() {
  const inputElement = document.getElementById("ms2-match-option-binsize-input") as HTMLInputElement;
  console.log(inputElement!.value);
  hideMatchPopup();

  let binsize: number = Number(inputElement!.value);
  let [names, identifiers, adducts, cids, mzs, cossim, matrices]: [string[], string[], string[], string[], string[], number[], string[]] = await invoke("match_msms_to_ui", {binsize:binsize});

  //console.log(identifiers, adducts, cids, mzs, cossim);
  let a: MSMSDatabase ={
    names: names,
    identifiers: identifiers,
    adducts: adducts,
    cids: cids,
    windows: [],
    tofs: [],
    mzs: mzs,
    cossim: cossim,
    matrices: matrices
  };
  updateMSMSResults(true, a);
  

}

export function match_msms() {

  MatchListener = () => {
    run_msms_match();
  }

  showMatchPopup();
  //run_msms_match();
}

/*
function addUserSpectra() {

	// Step 1: Select the parent container
	const parentContainer = document.querySelector("#ms2-results-body") as HTMLDivElement;

	if (parentContainer) {
		// Step 2: Create a new <div> element
		const template = document.createElement("template");


		template.innerHTML = generate_ms2_small_results_card("User Spectra", "_", "_", "_", 0.0, -1);
		const templateContent = template.content;


		// Step 4: Insert the new <div> element at the beginning of the parent container
		// You can use either appendChild or insertBefore, here's an example using insertBefore
		parentContainer.insertBefore(templateContent, parentContainer.firstChild);



    let a = document.getElementById("ms2-small-results-card-" + "-1") as HTMLDivElement;
		let b = a.querySelector("div.ms2-small-results-card-name") as HTMLDivElement;
		small_results_card_expanding(a, b, "USER1", -1, "0");
	}
}
*/