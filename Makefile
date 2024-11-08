setup:
	@chmod +x ./setup.sh
	./setup.sh

run:
	python main.py

clean:
	rm -rf ./coo *.jpg *.out *.png