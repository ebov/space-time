authors=[
{fn:"Gytis",ln:"Dudas",em:"gdudas@fredhutch.org",af:"Institute of Evolutionary Biology, University of Edinburgh, King's Buildings, Edinburgh, EH9 3FL, UK"},
//.
//. fill in authors here
//.
{fn:"Andrew",ln:"Rambaut",em:"a.rambaut@ed.ac.uk",af:"Institute of Evolutionary Biology, University of Edinburgh, King's Buildings, Edinburgh, EH9 3FL, UK"}
];
inputs=document.getElementsByTagName("input");
j=0;
for(i=0;i<inputs.length;i++){
	if (inputs[i].id=="authorEmail"+(j+1)) inputs[i].value=authors[j].em; 
	if (inputs[i].id=="authorFname"+(j+1)) inputs[i].value=authors[j].fn; 
	if (inputs[i].id=="authorLname"+(j+1)) inputs[i].value=authors[j].ln; 
	if (inputs[i].id=="authorAffil"+(j+1)) {
		inputs[i].value=authors[j].af;
		j++;
	}
}


