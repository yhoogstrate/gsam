
//___________________________________________Correlation of CD31 and PDGFRb___________________________________________//
// Markers:   COL1A1:Collagen, CD31:Endothelial Cells, PDGFR-B: Pericytes, GFAP:tumor cells and astrocytes

//rename image, split channels and close windows that are not needed
rename("sample");
run("Stack to Images");

selectWindow("FITC");
close();

selectWindow("Gold");
rename("Endothelial");

selectWindow("Cy5");
rename("Pericytes");

selectWindow("Red610");
close();

selectWindow("DAPI");
close();


//___________________________________________start of analysis___________________________________________//

//Measure average signal intensity of the pericyte marker
selectWindow("Pericytes");
run("Set Measurements...", "area mean standard modal min integrated display redirect=None decimal=3");
run("Measure");
Table.rename("Results", "Results-Mode-Pericytes");

//Measure average signal intensity of endothelial marker
selectWindow("Endothelial");
run("Set Measurements...", "area mean standard modal min integrated display redirect=None decimal=3");
run("Measure");
Table.rename("Results", "Results-Mode-Endothelial");
