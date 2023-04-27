//import { appWindow } from '@tauri-apps/api/window'
import { invoke } from '@tauri-apps/api';
import { open } from '@tauri-apps/api/dialog';
//import { fs } from "@tauri-apps/api/tauri";



window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("ms1-sidebar-open-file-button")!.addEventListener("click", async () => get_csv())
})


async function get_csv() {
    const result = await open({
        directory: false,
        multiple: false,
        filters: [{name: 'CSV Files', extensions:['csv']}],
    });
    
    if (result && result.length > 0) {
        const filePath = result;
        // Do something with the file path...

        let parsed_csv: string[] = await invoke("parse_ms1_csv", {path: filePath});

        console.log(filePath);
        
    }
  
}






