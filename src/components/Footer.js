import React from 'react';
import { Link } from 'react-router-dom';
import './Footer.css'; // Import the CSS file for styling

const Footer = () => (
  <footer className="footer">
    <div className="footer-content">
    <div className="footer-links">
        <h4>DrugForge</h4>
        <ul>
          <li><Link to="/">Home</Link></li>
          <li><Link to="/target-identification">Target Identification</Link></li>
          <li><Link to="/drug-properties-prediction">Drug Properties Prediction</Link></li>
          <li><Link to="/virtual-screening">Virtual Screening</Link></li>
          <li><Link to="/molecular-docking">Molecular Docking</Link></li>
        </ul>
      </div>
      <div className="footer-links">
        <h4>Services</h4>
        <ul>
          <li><Link to="/">Home</Link></li>
          <li><Link to="/target-identification">Target Identification</Link></li>
          <li><Link to="/drug-properties-prediction">Drug Properties Prediction</Link></li>
          <li><Link to="/virtual-screening">Virtual Screening</Link></li>
          <li><Link to="/molecular-docking">Molecular Docking</Link></li>
        </ul>
      </div>
      <div className="footer-contact">
        <h4>Contact Us</h4>
        <p>Email: info@drugforge.com</p>
        <p>Phone: +91 6382143070</p>
      </div>
      <div className="footer-social">
        <h4>Follow Us</h4>
        <ul className="social-links">
          {/* Use valid hrefs for external links */}
          <li><a href="https://www.facebook.com" target="_blank" rel="noopener noreferrer" aria-label="Facebook">Facebook</a></li>
          <li><a href="https://www.twitter.com" target="_blank" rel="noopener noreferrer" aria-label="Twitter">Twitter</a></li>
          <li><a href="https://www.linkedin.com" target="_blank" rel="noopener noreferrer" aria-label="LinkedIn">LinkedIn</a></li>
          <li><a href="https://www.instagram.com" target="_blank" rel="noopener noreferrer" aria-label="Instagram">Instagram</a></li>
        </ul>
      </div>
    </div>
    <div className="footer-bottom">
      <p>&copy; {new Date().getFullYear()} DrugForge. All rights reserved.</p>
    </div>
  </footer>
);

export default Footer;