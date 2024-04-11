from flask import Flask
import orbit_app as oapp

app = Flask(__name__)

@app.route("/")
def index():
    #m_E = oapp.m_Earth
    return "Hello world!"

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080, debug=True)