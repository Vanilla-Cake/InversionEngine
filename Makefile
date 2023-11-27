


all:
	cd src/Core/ && $(MAKE)
	cd src/Migration/GaussianBeam && $(MAKE)
	cd src/Tools && $(MAKE)
	cd src/Forward && $(MAKE)

install:
	cd src/Core/ && $(MAKE) install
	cd src/Migration/GaussianBeam && $(MAKE) install
	cd src/Tools && $(MAKE) install
	cd src/Forward && $(MAKE) install

clean:
	cd src/Core/ && $(MAKE) clean
	cd src/Migration/GaussianBeam && $(MAKE) clean
	cd src/Tools && $(MAKE) clean
	cd src/Forward && $(MAKE) clean
	


