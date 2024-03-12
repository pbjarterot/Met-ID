import { appWindow } from '@tauri-apps/api/window'
import './ms2/ms2_main';
import './ms1/ms1_main';
import { get_ctrl_v_data } from './ms1/ms1_popup';
import { invoke } from '@tauri-apps/api/tauri'
import { new_tgt_matrix } from './dropdown';

window.addEventListener("DOMContentLoaded", () => {
	invoke('close_splashscreen')
	document.getElementById('titlebar-minimize')!.addEventListener('click', () => appWindow.minimize());
	document.getElementById('titlebar-maximize')!.addEventListener('click', () => appWindow.toggleMaximize());
	document.getElementById('titlebar-close')!.addEventListener('click', () => appWindow.close());

	document.getElementById("tab-2")?.click();
	//document.getElementById("ms1-sidebar-add-metabolite")?.click();
	new_tgt_matrix()
});


document.addEventListener("keydown", async function(e) {
	if (e.ctrlKey && e.key === 'v') {
		const radioElements = document.querySelectorAll('input[type=radio]');
		// Initialize a variable to hold the value of the currently selected radio button
		let selectedValues: boolean[] | null = [];

		// Loop through and find which one is selected
		radioElements.forEach((radio: Element) => {
			const radioInput = radio as HTMLInputElement; // Typecast to HTMLInputElement
			selectedValues!.push(radioInput.checked); // Or radioInput.id if you want the ID of the selected element

		});

		// Now, `selectedValue` holds the value of the currently selected radio button

		console.log(`Currently selected tab is: ${selectedValues[1]}`);

		if (selectedValues[1] === true) {
			try {
				const clipboardText = await navigator.clipboard.readText();
				get_ctrl_v_data(clipboardText);
				console.log('Pasted content:', clipboardText);
			  // Now you can send clipboardText to the Tauri backend if needed
			} catch (err) {
			  console.error('Could not read from clipboard:', err);
			}
		}
		
	  }
});




