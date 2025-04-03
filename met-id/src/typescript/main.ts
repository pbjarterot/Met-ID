import { getCurrentWebviewWindow } from '@tauri-apps/api/webviewWindow'
import './ms2/ms2_main';
import './ms1/ms1_main';
import './database_tab/db_main.ts';
import { invoke } from '@tauri-apps/api/core'
import { listen } from "@tauri-apps/api/event";
import { new_tgt_matrix } from './dropdown';
const appWindow = getCurrentWebviewWindow()






window.addEventListener("DOMContentLoaded", () => {
	invoke('close_splashscreen')
	document.getElementById('titlebar-minimize')!.addEventListener('click', () => appWindow.minimize());
	document.getElementById('titlebar-maximize')!.addEventListener('click', () => appWindow.toggleMaximize());
	document.getElementById('titlebar-close')!.addEventListener('click', () => appWindow.close());

	document.getElementById("tab-2")?.click();
	new_tgt_matrix()
	
	listen<string>("ask-frontend-bool", async (event) => {
		const id = event.payload;
		const userResponse = window.confirm("Would you like to update Met-ID?");
	
		await invoke("frontend_bool_response", {
		id,
		value: userResponse,
		});
	});
});
