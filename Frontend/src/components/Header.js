// Header.js
import React from 'react';
import { Link, useLocation } from 'react-router-dom';
import './Header.css';

const Header = () => {
  const location = useLocation();

  return (
    <header className="navbar">
      <div className="navbar-logo">
        <Link to="/" className="logo">DrugForge</Link>
      </div>
      <nav className="navbar-links">
        <ul>
          <li>
            <Link to="/services" className={location.pathname === '/services' ? 'active' : ''}>
              Services
            </Link>
          </li>
          <li>
            <Link to="/blog" className={location.pathname === '/blog' ? 'active' : ''}>
              Blogs
            </Link>
          </li>
          <li>
            <Link to="/features" className={location.pathname === '/features' ? 'active' : ''}>
              Features
            </Link>
          </li>
          <li>
            <Link to="/contact" className={location.pathname === '/contact' ? 'active' : ''}>
              Contact Us
            </Link>
          </li>
          <li>
            <Link to="/pricing" className={location.pathname === '/pricing' ? 'active' : ''}>
              Pricing
            </Link>
          </li>
        </ul>
      </nav>
      <div className="navbar-auth">
        <Link to="/signin" className="sign-in-button">Sign In</Link>
        <Link to="/profile" className="profile-button">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="profile-icon">
            <path fillRule="evenodd" d="M18.685 19.097A9.723 9.723 0 0021.75 12c0-5.385-4.365-9.75-9.75-9.75S2.25 6.615 2.25 12a9.723 9.723 0 003.065 7.097A9.716 9.716 0 0012 21.75a9.716 9.716 0 006.685-2.653zm-12.54-1.285A7.486 7.486 0 0112 15a7.486 7.486 0 015.855 2.812A8.224 8.224 0 0112 20.25a8.224 8.224 0 01-5.855-2.438zM15.75 9a3.75 3.75 0 11-7.5 0 3.75 3.75 0 017.5 0z" clipRule="evenodd" />
          </svg>
          Profile
        </Link>
        <Link to="/dashboard" className="profile-button">
          Dashboard
        </Link>
      </div>
    </header>
  );
};

export default Header;