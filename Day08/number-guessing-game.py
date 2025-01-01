import tkinter as tk
from tkinter import messagebox
import random

class NumberGuessingGame:
    def __init__(self, root):
        self.root = root
        self.root.title("Number Guessing Game")
        self.root.geometry("400x300")
        
        self.number = random.randint(1, 20)
        self.attempts = 0
        
        # Create GUI elements
        self.label = tk.Label(root, text="Enter a number between 1 and 20!", font=("Arial", 12))
        self.label.pack(pady=20)
        
        self.attempts_label = tk.Label(root, text="Attempts: 0", font=("Arial", 10))
        self.attempts_label.pack()
        
        self.entry = tk.Entry(root, font=("Arial", 14))
        self.entry.pack(pady=20)
        
        self.submit_button = tk.Button(root, text="Submit Guess", command=self.handle_guess)
        self.submit_button.pack(pady=5)
        
        self.show_button = tk.Button(root, text="Show Number", command=self.show_number)
        self.show_button.pack(pady=5)
        
        self.new_game_button = tk.Button(root, text="New Game", command=self.new_game)
        self.new_game_button.pack(pady=5)
        
        # Bind Enter key to submit
        self.entry.bind('<Return>', lambda event: self.handle_guess())
        
    def handle_guess(self):
        guess = self.entry.get().strip().lower()
        
        if guess.isdigit():
            number_guess = int(guess)
            self.attempts += 1
            self.attempts_label.config(text=f"Attempts: {self.attempts}")
            
            if number_guess > self.number:
                message = "Your number is too big. Enter another."
            elif number_guess < self.number:
                message = "Your number is too small. Enter another."
            else:
                message = f"Exactly! The correct number is {self.number}.\nYou guessed it in {self.attempts} attempts."
                messagebox.showinfo("Congratulations!", message)
                self.disable_game()
            
            self.label.config(text=message)
        else:
            messagebox.showerror("Error", "Please enter a valid number between 1 and 20!")
        
        self.entry.delete(0, tk.END)
        self.entry.focus()
    
    def show_number(self):
        messagebox.showinfo("Number Revealed", f"The hidden number is {self.number}")
    
    def new_game(self):
        self.number = random.randint(1, 20)
        self.attempts = 0
        self.attempts_label.config(text="Attempts: 0")
        self.label.config(text="Enter a number between 1 and 20!")
        self.enable_game()
        self.entry.delete(0, tk.END)
        self.entry.focus()
    
    def disable_game(self):
        self.entry.config(state='disabled')
        self.submit_button.config(state='disabled')
        self.show_button.config(state='disabled')
    
    def enable_game(self):
        self.entry.config(state='normal')
        self.submit_button.config(state='normal')
        self.show_button.config(state='normal')

if __name__ == '__main__':
    root = tk.Tk()
    game = NumberGuessingGame(root)
    root.mainloop()
