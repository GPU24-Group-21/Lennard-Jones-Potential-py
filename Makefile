setup:
	@chmod +x ./setup.sh
	./setup.sh

run:
	@chmod +x ./run.sh
	./run.sh $(ARGS)

clean:
	rm -rf ./coo *.jpg *.out *.png