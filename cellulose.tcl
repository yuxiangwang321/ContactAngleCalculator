#  cellulose.tcl, execute multiple estimations  

source cac-cg.tcl
                             
cd bcwt
	cac bcwt.psf bcwt.dcd xz 0 49 1 {index 51273 to 62051} 3 0.0324 0.1 3 0.5 0.5
cd ..

cd L0 
	cac L0.psf L0.dcd xz 0 49 1 {index 53298 to 64076} 3 0.0324 0.1 3 0.5 0.5
cd .. 

cd L4
	cac L4.psf L4.dcd xz 0 49 1 {index 58158 to 68936} 3 0.0324 0.1 3 0.5 0.5
cd .. 

cd L8
	cac L8.psf L8.dcd xz 0 49 1 {index 63018 to 73796} 3 0.0324 0.1 3 0.5 0.5
cd .. 

puts “ALL DONE!”
