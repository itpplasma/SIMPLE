VMEC=xvmec

wout.nc: input.purely_toroidal_field
	${VMEC} purely_toroidal_field
	mv wout_purely_toroidal_field.nc wout.nc

input.purely_toroidal_field:
	wget https://raw.githubusercontent.com/hiddenSymmetries/simsopt/refs/heads/master/tests/test_files/input.purely_toroidal_field 

clean:
	rm -f dcon_purely_toroidal_field.txt input.purely_toroidal_field \
	log.purely_toroidal_field mercier.purely_toroidal_field \
	parvmecinfo.txt threed1.purely_toroidal_field timings.txt
