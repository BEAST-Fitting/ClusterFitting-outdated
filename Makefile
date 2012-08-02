.SUFFIXES: .fits.gz

all:
	make data

data: data/foo.fits

%.fits.gz:
	curl -o data/$*.fits "http://karl.com/$*.fits"
