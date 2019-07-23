import os
import pandas as pd

from flask import (Flask, 
        request, 
        render_template, 
        url_for, 
        redirect, 
        jsonify, session)

#Formatting the file name
from werkzeug.utils import secure_filename

# Exceptions
from werkzeug.exceptions import RequestEntityTooLarge

# Import csrf token
from flask_wtf.csrf import CSRFProtect 

from env import *

#Import Controller
from controllers.fileUploadController import *



#-----------------------------------------------------------
#Initialization section
#-----------------------------------------------------------
app = Flask(__name__)
csrf = CSRFProtect(app)




#-----------------------------------------------------------
# Config Constants
#-----------------------------------------------------------
app.config['SECRET_KEY'] = SECRET_KEY
app.config['FILE_UPLOAD_PATH'] = FILE_UPLOAD_PATH
app.config['ALLOWED_FILES'] = ALLOWED_FILES

app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH

# app.config['MAX_SIZE_FILE'] = MAX_SIZE_FILE



#-----------------------------------------------------------
# Helper function
#-----------------------------------------------------------

#Check if the file name is correct
def allowed_extensions(filename):
    if "." not in filename:
        return False
    #Get the extension of the file
    ext = filename.rsplit(".", 1)[-1]
    if ext.upper() in app.config['ALLOWED_FILES']:
        return True
    
    return False

#-----------------------------------------------------------
# Routes 
#-----------------------------------------------------------
#######################
#Render the index page
#######################
@app.route('/')
def index():    
    return render_template("index.html", title="Welcome")

#######################################################
# Upload file, parse it and return the list of samples
########################################################

@app.route('/file_upload', methods=['POST'])
def upload_file():
# See the following post for step by step
# https://www.shanelynn.ie/asynchronous-updates-to-a-webpage-with-flask-and-socket-io/
    if request.method == "POST":
        errors = []
        fichier = ""
        if request.files:

            try :
                fichier = request.files['input_file']
            except RequestEntityTooLarge as e:
                return jsonify(status="error", messages= ["file too large"])

            # Check in the file has a name
            if fichier.filename == "" or fichier.filename is None:
                errors.append("The file has no name.")

            #Check if the extension of the file is the one expected
            if not allowed_extensions(fichier.filename):
                errors.append("The file extension you provided is not supported.")            

            if len(errors) > 0 :                
                return jsonify(status="error", messages= errors)

            else:
                #Create a secure fileName
                secure_fname = secure_filename(fichier.filename)
                path_to_file = os.path.join(app.config['FILE_UPLOAD_PATH'], secure_fname)               

                # Save the file                
                fichier.save(path_to_file) 

                # Validating the content of the file
                errors = validate_uploaded_file(path_to_file)
                if len(errors) > 0 :
                    return jsonify(status="error", messages=errors)           


                # Save the path of the file to session since we can't store the full dataset is session
                session["data_set_path"] = path_to_file  
               
                
                #Return the sample names               
                # data, taxonomy = return_sample_names_and_taxonomy(session["data_set_path"])
                data = return_sample_names(session["data_set_path"])

                return jsonify(status="success", 
                                messages=["Successfully uploaded"],                                
                                samples=data,
                                # taxonomy_detailed= taxonomy
                                )
        else:
            return "Invalid input provides"
    else:
        return "An unknown error occured"


##################################
# Return taxa abundance by sample
##################################
@app.route('/taxonomy_by_sample_name/<sample>')    
def taxonomy_by_sample_name(sample):
    '''
    Return the taxonomy: taxa with non zero counts based on the selected sample
    zero count taxa are not necessary especially for plotting, they will just add elements to the legend
    '''
    taxonomy = return_taxonomy_by_sample_name(session["data_set_path"], sample)    
    return jsonify(status= "success", taxonomy = taxonomy)


##################################
# Return taxa abundance by sample
##################################
@app.route('/abundance_by_sample/<sample>')
def abundance_by_sample(sample):
    abundance = return_abundance_by_taxonomy_level(session["data_set_path"], sample, 2)
    print(abundance)

    return jsonify(status="success", 
                        messages=["Sample fetched successfully"],                         
                        abundance= abundance)



##################################
# Return taxa abundance by sample
##################################
@app.route('/taxon_abundance_by_sample/<sample>/<level>/<taxon>')
def abundance_taxon_by_sample(sample, level, taxon):
    abundance = return_abundance_by_taxon_by_sample(session["data_set_path"], sample, level, taxon)
    # print(abundance)

    return jsonify(status="success", 
                        messages=["Sample fetched successfully"],                         
                        abundance= abundance)





#------------------------------------
# Handeling errors
#------------------------------------

@app.errorhandler(413)
@app.errorhandler(RequestEntityTooLarge)
def app_handle_413(e):
    return 'File Too Large'


if __name__ == "__main__":
    app.run(debug=True)

