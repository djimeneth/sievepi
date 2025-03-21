all:		sievepi

sievepi:	sievepi.cc
		g++ -Wall -O3 -fomit-frame-pointer -o sievepi sievepi.cc
