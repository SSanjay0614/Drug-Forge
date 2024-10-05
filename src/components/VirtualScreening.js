import React, { useState } from 'react';

function VirtualScreening() {
  const [molecule, setMolecule] = useState('');
  const [screeningResults, setScreeningResults] = useState([]);

  const handleScreenMolecule = (e) => {
    e.preventDefault();
    
    // Simulate screening results
    const results = [
      { name: molecule, score: Math.random() * 100 },
      { name: molecule + ' Variant A', score: Math.random() * 100 },
      { name: molecule + ' Variant B', score: Math.random() * 100 },
    ];
    
    setScreeningResults(results);
  };

  return (
    <div>
      <h1>Virtual Screening</h1>
      <form onSubmit={handleScreenMolecule}>
        <label>
          Molecule Name:
          <input
            type="text"
            value={molecule}
            onChange={(e) => setMolecule(e.target.value)}
            required
          />
        </label>
        <button type="submit">Screen Molecule</button>
      </form>

      {screeningResults.length > 0 && (
        <div>
          <h2>Screening Results</h2>
          <ul>
            {screeningResults.map((result, index) => (
              <li key={index}>
                {result.name}: Score {result.score.toFixed(2)}
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}

export default VirtualScreening;