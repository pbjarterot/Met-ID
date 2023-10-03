import { All_spectra_database } from './ms2_main';

export function generate_result_card(spectra: All_spectra_database, i: number) {
const template = `
		<div class="ms2-results-card">
			<div class="ms2-results-card-textbox">
				<div class="ms2-results-card-name">
				<label class="ms2-results-card-name">${spectra[0][i]}</label>
				</div>
				<div class="ms2-results-card-dash">
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>CID</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>${spectra[2][i]}</p>
					</div>
				</div>
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>Matrix</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>FMP-10</p>
					</div>
				</div>
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>Adduct</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>${spectra[3][i]}</p>
					</div>
				</div>
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>Parent ion m/z</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>400.000000</p>
					</div>
				</div>
				</div>
				<div class="molecule" id="molecule-canvas">
				<img data-smiles=${spectra[1][i]}
						data-smiles-options=" {'width': 300, 'height': 250 }"  
						data-smiles-theme='dark' />
				
				<script>
					SmiDrawer.apply({"themes": dark});
				</script>

				</div>
			</div>
			<div class="ms2-results-card-spectrogram-div" id="ms2-results-card-spectrogram-div${i as unknown as string}"></div>

			</div>
            `
    return template
}