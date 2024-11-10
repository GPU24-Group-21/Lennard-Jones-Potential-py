if [ ! -d "venv" ]; then
    # create virtual environment
    python3 -m venv venv
    # activate virtual environment
    source venv/bin/activate
    # install requirements
    pip install -r requirements.txt
else
    # activate virtual environment
    source venv/bin/activate
fi

# install requirements
pip install -r requirements.txt