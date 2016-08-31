hd = $(HOME)/cosmo/lib
LIB = -lm -fopenmp -L${hd} -lcutil 


#CC = gcc
CC = clang
CFLAGS = -O2


OBJS10 = cen_sat_decomp.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o
cen_sat_decomp:	$(OBJS10)
	$(CC) -o $@ $(OBJS10) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS11 = cen_sat_test.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o
cen_sat_test:	$(OBJS11)
	$(CC) -o $@ $(OBJS11) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS12 = galaxy_environments.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o
galaxy_environments:	$(OBJS12)
	$(CC) -o $@ $(OBJS12) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS13=	construct_SHAM_hod.o  meshlink2.o nbrsfind2.o i3tensor_2.o
construct_SHAM_hod:	$(OBJS13)
	$(CC) -o $@ $(OBJS13) $(LIB)
	cp $@ $(HOME)/exec/

OBJS14= 2Dbinning.o
2Dbinning:	$(OBJS14)
	$(CC) -o $@ $(OBJS14) $(LIB)
	cp $@ $(HOME)/exec/


OBJS15 = clf_group_finder.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder:	$(OBJS15)
	$(CC) -o $@ $(OBJS15) $(LIB)
	cp -f $@ $(HOME)/exec/$@


OBJS16 = clf_group_finder_mocks.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o ran1.o gasdev.o 
clf_group_finder_mocks:	$(OBJS16)
	$(CC) -o $@ $(OBJS16) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS16a = clf_group_finder_mocks_old1.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o ran1.o gasdev.o 
clf_group_finder_mocks_old1:	$(OBJS16a)
	$(CC) -o $@ $(OBJS16a) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS17 = identify_hostmass.o 
identify_hostmass:	$(OBJS17)
	$(CC) -o $@ $(OBJS17) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS18 = galaxy_environments_random.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o
galaxy_environments_random:	$(OBJS18)
	$(CC) -o $@ $(OBJS18) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS181 = galaxy_environments_random2.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o
galaxy_environments_random2:	$(OBJS181)
	$(CC) -o $@ $(OBJS181) $(LIB)
	cp -f $@ $(HOME)/exec/$@


OBJS19 = clf_group_finder_mstellar.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder_mstellar:	$(OBJS19)
	$(CC) -o $@ $(OBJS19) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS20 = nearest_group.o midpnt.o qromo.o polint.o
nearest_group:	$(OBJS20)
	$(CC) -o $@ $(OBJS20) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS201 = conformity_groups.o midpnt.o qromo.o polint.o sort2.o
conformity_groups:	$(OBJS201)
	$(CC) -o $@ $(OBJS201) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS203 = conformity_groups_kauffmann.o midpnt.o qromo.o polint.o sort2.o sort.c
conformity_groups_kauffmann:	$(OBJS203)
	$(CC) -o $@ $(OBJS203) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS204 = conformity_groups_sfr.o midpnt.o qromo.o polint.o sort2.o sort.c
conformity_groups_sfr:	$(OBJS204)
	$(CC) -o $@ $(OBJS204) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS202 = conformity_groups_mock.o midpnt.o qromo.o polint.o sort2.o
conformity_groups_mock:	$(OBJS202)
	$(CC) -o $@ $(OBJS202) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS21 = galaxy_environments_nonlocal.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o
galaxy_environments_nonlocal:	$(OBJS21)
	$(CC) -o $@ $(OBJS21) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS22 = clf_group_finder_primus.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o
clf_group_finder_primus:	$(OBJS22)
	$(CC) -o $@ $(OBJS22) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS23 = mstellar_subset.o
mstellar_subset:	$(OBJS23)
	$(CC) -o $@ $(OBJS23) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS24 = galaxy_environments_mstellar.o spline.o splint.o qromo.o midpnt.o \
	polint.o sort2.o sham.o zbrent.o
galaxy_environments_mstellar:	$(OBJS24)
	$(CC) -o $@ $(OBJS24) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS25 = clf_group_finder_mocks_mstellar.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o ran1.o gasdev.o 
clf_group_finder_mocks_mstellar:	$(OBJS25)
	$(CC) -o $@ $(OBJS25) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS26 = clf_group_finder_bolshoi.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder_bolshoi:	$(OBJS26)
	$(CC) -o $@ $(OBJS26) $(LIB)
	cp -f $@ $(HOME)/exec/$@


OBJS27 = clf_group_finder_bolshoi_radec.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder_bolshoi_radec:	$(OBJS27)
	$(CC) -o $@ $(OBJS27) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS271 = clf_group_finder_bolshoi_radec_omp.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder_bolshoi_radec_omp:	$(OBJS271)
	$(CC) -o $@ $(OBJS271) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS28 = nearest_group_ms.o midpnt.o qromo.o polint.o
nearest_group_ms:	$(OBJS28)
	$(CC) -o $@ $(OBJS28) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS29 = clf_group_finder_nsa.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder_nsa:	$(OBJS29)
	$(CC) -o $@ $(OBJS29) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS30 = clf_group_finder_nsa_vmax.o spline.o splint.o qromo.o midpnt.o polint.o sort2.o sham.o zbrent.o \
	trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
clf_group_finder_nsa_vmax:	$(OBJS30)
	$(CC) -o $@ $(OBJS30) $(LIB)
	cp -f $@ $(HOME)/exec/$@



clean:
	rm -f *.o
