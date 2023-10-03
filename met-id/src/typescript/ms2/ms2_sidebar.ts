import { open } from '@tauri-apps/api/dialog';

let result: string | string[] | null = "";

export async function get_msms_from_mzml() {
    result = await open({directory: false, multiple: false, filters: [{name: 'mzML', extensions:['mzML']}]});

    if (result && result.length > 0) {


    }
}