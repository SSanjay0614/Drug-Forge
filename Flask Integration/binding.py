from flask import Flask, request, jsonify
import pickle
from DeepPurpose import utils
from DeepPurpose import DTI as models
from flask_cors import CORS  # For handling Cross-Origin Resource Sharing

# Initialize Flask application
app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})  # Enable CORS for all routes

# Load the pre-trained model
with open('binding_model.pkl', 'rb') as f:
    model = pickle.load(f)

# Route to check API status
@app.route('/', methods=['GET'])
def get_data():
    return jsonify({"message": "Drug-Target Interaction Model API is Running"})

# Define the route for making predictions
@app.route('/predict', methods=['POST'])
def predict():
    try:
        # Get the drug and target from the POST request
        data = request.get_json()
        drug = data.get('drug')
        target = data.get('target')

        if not drug or not target:
            return jsonify({'error': 'Drug or target is missing'}), 400

        # Create input for the model
        X_drug = [drug]
        X_target = [target]
        y = [0]  # Dummy value, not used in prediction
        drug_encoding, target_encoding = 'MPNN', 'CNN'
        X_pred = utils.data_process(X_drug, X_target, y, 
                                    drug_encoding, target_encoding, 
                                    split_method='no_split')

        # Model prediction
        y_pred = model.predict(X_pred)

        # Return the prediction as a JSON response
        return jsonify({'predicted_score': y_pred[0]})

    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Run the Flask application
if __name__ == '__main__':
    app.run(debug=True, port=5004)
