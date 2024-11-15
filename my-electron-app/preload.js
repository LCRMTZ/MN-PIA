const { contextBridge } = require('electron');

contextBridge.exposeInMainWorld('electronAPI', {
    addEventListenerToButton: (callback) => {
        document.addEventListener('DOMContentLoaded', () => {
            const button = document.getElementById('compare-button');
            if (button) {
                button.addEventListener('click', callback);
            }
        });
    }
});
