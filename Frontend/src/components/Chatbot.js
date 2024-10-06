import React, { useState } from 'react';
import './Chatbot.css'; // Import CSS file for styling
import axios from 'axios';
import ReactMarkdown from 'react-markdown';

const Chatbot = () => {
    const [isOpen, setIsOpen] = useState(false);
    const [message, setMessage] = useState('');
    const [messages, setMessages] = useState([
        { text: 'Welcome to DrugForge! How can I assist you today?', sender: 'bot' }
    ]);
    const [generatingAnswer, setGeneratingAnswer] = useState(false);

    const toggleChat = () => {
        setIsOpen(!isOpen);
    };

    const closeChat = () => {
        setIsOpen(false);
    };

    const handleSendMessage = async () => {
        if (message.trim() !== '') {
            setMessages([...messages, { text: message, sender: 'user' }]);
            setMessage('');
            setGeneratingAnswer(true);
            try {
                const response = await axios({
                    url: `https://generativelanguage.googleapis.com/v1beta/models/gemini-pro:generateContent?key=${
                        import.meta.env.VITE_API_GENERATIVE_LANGUAGE_CLIENT
                    }`,
                    method: "post",
                    data: {
                        contents: [{ parts: [{ text: message }] }],
                    },
                });

                setMessages([...messages, { text: response["data"]["candidates"][0]["content"]["parts"][0]["text"], sender: 'bot' }]);
            } catch (error) {
                console.log(error);
                setMessages([...messages, { text: "Sorry - Something went wrong. Please try again!", sender: 'bot' }]);
            }

            setGeneratingAnswer(false);
        }
    };

    return (
        <div className="chatbot-container">
            {isOpen ? (
                <div className="chatbot-window">
                <div className="chatbot-header">
                    <h3>Chat with Us</h3>
                    <button className="close-button" onClick={closeChat}>✖</button>
                </div>
                <div className="chatbot-body">
                    <div className="messages">
                        {messages.map((message, index) => (
                            <div key={index} className={`message ${message.sender === 'bot' ? 'bot-message' : 'user-message'}`}>
                                <div className="message-header">
                                    <span className="message-sender">{message.sender === 'bot' ? 'Bot' : 'You'}</span>
                                    <span className="message-timestamp">{new Date().toLocaleTimeString()}</span>
                                </div>
                                <div className="message-content">
                                    <ReactMarkdown className="p-4">{message.text}</ReactMarkdown>
                                </div>
                            </div>
                        ))}
                    </div>
                </div>
                <div className="chatbot-input">
                    <input
                        type="text"
                        placeholder="Type your message..."
                        value={message}
                        onChange={(e) => setMessage(e.target.value)}
                        className="input-field"
                    />
                    <button onClick={handleSendMessage} disabled={generatingAnswer} className="send-button">
                        {generatingAnswer ? 'Generating answer...' : 'Send'}
                    </button>
                </div>
            </div>
            ) : (
                <div className="chatbot-icon">
                    <video width="50" height="50" loop autoPlay muted>
                        <source src="https://res.cloudinary.com/dvude7m7p/video/upload/v1728057294/DrugForge/kdoaxazwcl8tqx30irun.mp4" type="video/mp4" />
                        Your browser does not support the video tag.
                    </video>
                    <button onClick={toggleChat}> O </button>
                </div>
            )}
        </div>
    );
};

export default Chatbot;