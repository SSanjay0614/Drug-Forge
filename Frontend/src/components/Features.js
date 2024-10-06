import React from 'react';
import { Link } from 'react-router-dom';
import './Features.css'; // Import the corresponding CSS file

const features = [
  {
    title: 'Solubility Prediction',
    description: 'Predict the solubility of chemical compounds.',
    imageUrl: 'https://res.cloudinary.com/dvude7m7p/image/upload/v1727851104/DrugForge/ascbt6rplkp8ghczkvl9.jpg', // Replace with actual image path
    route: '/predictor',
  },
  {
    title: 'Blood-Brain Barrier Permeability',
    description: 'Evaluate the permeability of compounds across the blood-brain barrier.',
    imageUrl: 'https://res.cloudinary.com/dvude7m7p/image/upload/v1727851104/DrugForge/txupp5fotqaz0shdxfyv.jpg', // Replace with actual image path
    route: '/bbb-permeability', // Change this route as needed
  },
  {
    title: 'Target Identification and Activity Prediction',
    description: 'Identify potential drug targets and predict their activity.',
    imageUrl: 'https://res.cloudinary.com/dvude7m7p/image/upload/v1727851104/DrugForge/n7b3njatrknrdlh1dm8k.jpg', // Replace with actual image path
    route: '/target-identification', // Change this route as needed
  },
  {
    title: 'Molecular Docking',
    description: 'Analyze the interaction between molecules and targets.',
    imageUrl: 'https://res.cloudinary.com/dvude7m7p/image/upload/v1727851104/DrugForge/avuisc187v0k54bofhni.jpg', // Replace with actual image path
    route: '/molecular-docking', // Change this route as needed
  },
];

const Features = () => {
  return (
    <div className="features-container">
      <h2>Our Features</h2>
      <div className="features-grid">
        {features.map((feature, index) => (
          <Link key={index} to={feature.route} className="feature-card">
            <img src={feature.imageUrl} alt={feature.title} className="feature-image" />
            <h3>{feature.title}</h3>
            <p>{feature.description}</p>
          </Link>
        ))}
      </div>
    </div>
  );
};

export default Features;
