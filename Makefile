setup:
	@chmod +x ./setup.sh
	./setup.sh

run:
	source ./venv/bin/activate && python ./main.py

clean:
	rm -rf ./coo *.jpg *.out *.png