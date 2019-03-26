from flask import send_file, jsonify
from flask_mail import Message
import re

def downloadFile(files, request):
    try:
        if (request.values['fileName'] in files):
            return send_file(files[request.values['fileName']])
        else:
            raise(Exception)
    except Exception as e:
        return jsonify(error_message='Your file was not found, it may have reached a time-out limit. Please re-run NAPS.')

def emailFiles(files, request, mail, app):
    if files:
        emailAddress = request.form['emailAddress']
        if not re.fullmatch(r"[^@]+@[^@]+\.[^@]+", emailAddress):
            return jsonify(status='email_failed', message='Invalid email address entered.')

        msg = Message("NAPS Results",
                sender=app.config.get("MAIL_USERNAME"),
                recipients=[emailAddress])

        msg.body = "NAPS results"

        if 'results' in files:
            with app.open_resource(files['results']) as fp:
                msg.attach("results.txt", "text/plain", fp.read())
        if 'plot' in files:
            with app.open_resource(files['plot']) as fp:
                msg.attach("plot.html", "text/html", fp.read())
        mail.send(msg)
        return jsonify(status='ok')
    else:
        return jsonify(status='email_failed', message='Email sending failed, please refresh and retry.')