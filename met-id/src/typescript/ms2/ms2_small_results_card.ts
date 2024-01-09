import { draw } from "./ms2_draw_spectra";
import * as d3 from 'd3';



export function generate_ms2_small_results_card(name: string, _cid: string, adduct: string, parent_mz: string | number, cossim: number, matrix: string, index: number) {
    if (!cossim) {
        cossim = 0.0;
    } 
    console.log(parent_mz, matrix)
    let template = `<div class="ms2-small-results-card" id="ms2-small-results-card-${index}">
                        <div class="ms2-small-results-card-top">
                            <div class="ms2-small-results-card-textbox">
                                <div class="ms2-small-results-card-name">
                                    <label class="ms2-small-results-card-name">${name}</label>
                                </div>
                            </div>
                            <div class="ms2-small-results-card-metadata">
                                <div class="ms2-small-results-molecule-info">
                                    <div class="small-ms2-dash-box">
                                        <div class="small-ms2-dashbox-header">
                                            <p>Cosine</p>
                                        </div>
                                        <div class="small-ms2-dashbox-text">
                                            <p>${(cossim as number).toFixed(3)}</p>
                                        </div>
                                    </div>
                                    <div class="small-ms2-dash-box">
                                        <div class="small-ms2-dashbox-header">
                                            <p>Adduct</p>
                                        </div>
                                        <div class="small-ms2-dashbox-text">
                                            <p>${adduct}</p>
                                        </div>
                                    </div>
                                    <div class="small-ms2-dash-box">
                                        <div class="small-ms2-dashbox-header">
                                            <p>Matrix</p>
                                        </div>
                                        <div class="small-ms2-dashbox-text">
                                            <p>${matrix}</p>
                                        </div>
                                    </div>

                                    <div class="small-ms2-dash-box">
                                        <div class="small-ms2-dashbox-header">
                                            <p>MS1 m/z</p>
                                        </div>
                                        <div class="small-ms2-dashbox-text">
                                            <p>${parent_mz}</p>
                                        </div>
                                    </div>
                                </div> 
                            </div>
                        </div>
                    </div>`

    return template
}

export function small_results_card_expanding(expandableDiv: HTMLDivElement, listeningDiv: HTMLDivElement, identifier: string, index: number, adduct: string, cid: string) {
    listeningDiv.addEventListener('click', () => {
        // Check if the div is already expanded
        if (expandableDiv.classList.contains('open')) {
            // Collapse the div
            expandableDiv.classList.remove('open');
            // Remove additional divs
            const additionalDiv = document.getElementById('ms2-small-results-card-plot-' + index.toString());
            additionalDiv!.remove();
            //additionalDivs.forEach((div) => div.remove());
        } else {
            // Expand the div
            expandableDiv.classList.add('open');

            const newDiv = document.createElement('div');
            newDiv.className = 'ms2-small-results-card-plot';
            newDiv.id = 'ms2-small-results-card-plot-' + index.toString();

            expandableDiv.appendChild(newDiv);
            draw(index, identifier, adduct, cid);

            // Select the SVG container
            var svgContainer = d3.select("#ms2-small-results-card-plot-" + index.toString() as unknown as string);

            // Add event listener for mouseover
            svgContainer.on('mouseover', function () {
            window.addEventListener('wheel', preventScrolling, { passive: false });
            });

            // Add event listener for mouseout
            svgContainer.on('mouseout', function () {
            window.removeEventListener('wheel', preventScrolling);
            });

            // Function to prevent scrolling
            function preventScrolling(e) {
            e.preventDefault();
            }
        }
    });
}

