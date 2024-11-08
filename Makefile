setup:
	@chmod +x ./setup.sh
	./setup.sh

run:
	@chmod +x ./run.sh
	./run.sh

clean:
	rm -rf ./coo *.jpg *.out *.png