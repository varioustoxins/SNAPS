from flask import jsonify
from flask_mail import Message
import re

def emailFiles(request, mail, app):
    files = request.values
    if files:
        emailAddress = request.form['emailAddress']
        if not re.fullmatch(r"[^@]+@[^@]+\.[^@]+", emailAddress):
            return jsonify(status='email_failed', message='Invalid email address entered.')

        msg = Message("NAPS Results",
                sender=app.config.get("MAIL_USERNAME"),
                recipients=[emailAddress])

        msg.body = "NAPS results"

        if 'results' in files:
            msg.attach("results.txt", "text/plain", files['results'])
        if 'plot' in files:
            msg.attach("plot.html", "text/html", files['plot'])
        mail.send(msg)
        return jsonify(status='ok')
    else:
        return jsonify(status='email_failed', message='Email sending failed, please refresh and retry.')