import React from 'react';
import ReactDOM from 'react-dom/client'; // Use this import for React 18+
import App from './App';
import { BrowserRouter } from 'react-router-dom';

const root = ReactDOM.createRoot(document.getElementById('root')); // Create a root for rendering
root.render(
  <BrowserRouter>
    <App />
  </BrowserRouter>
);