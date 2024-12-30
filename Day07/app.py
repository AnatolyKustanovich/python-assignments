from flask import Flask, render_template, request, session, jsonify
import random

app = Flask(__name__)
app.secret_key = 'your-secret-key-here'

@app.route('/', methods=['GET', 'POST'])
def game():
    if 'number' not in session:
        initialize_game()
    
    return render_template('game.html', 
                         message=session.get('message', "Enter a number between 1 and 20!"),
                         attempts=session.get('attempts', 0),
                         game_over=session.get('game_over', False))

@app.route('/guess', methods=['POST'])
def handle_guess():
    if session.get('game_over', False):
        return jsonify({'message': 'Game is over. Start a new game!', 'game_over': True})
    
    guess = request.form.get('guess', '').strip().lower()
    
    if guess == 's':
        return jsonify({
            'message': f"The hidden number is {session['number']}.",
            'game_over': False
        })
    
    if guess.isdigit():
        number_guess = int(guess)
        session['attempts'] += 1
        
        if number_guess > session['number']:
            message = "Your number is too big. Enter another."
        elif number_guess < session['number']:
            message = "Your number is too small. Enter another."
        else:
            message = f"Exactly! The correct number is {session['number']}. "\
                     f"You guessed it in {session['attempts']} attempts."
            session['game_over'] = True
            
        return jsonify({
            'message': message,
            'attempts': session['attempts'],
            'game_over': session['game_over']
        })
    
    return jsonify({
        'message': "Invalid input! Please enter a number between 1 and 20, or 's' to show the number.",
        'game_over': False
    })

@app.route('/new-game', methods=['POST'])
def new_game():
    initialize_game()
    return jsonify({
        'message': "New game started! Enter a number between 1 and 20!",
        'attempts': 0,
        'game_over': False
    })

def initialize_game():
    session['number'] = random.randint(1, 20)
    session['attempts'] = 0
    session['game_over'] = False
    session['message'] = "Enter a number between 1 and 20!"

if __name__ == '__main__':
    app.run(debug=True)