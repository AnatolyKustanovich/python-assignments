
import random

def main():
    while True:
        play_game()
        if not ask_play_again():
            break
    print("Thanks for playing! Goodbye!")

def play_game():
    number_guess = random.randint(1, 20)
    guesses = 0
    print("Try to guess number between 1 and 20."  )
    while True:
        number_user = input()
        if number_user == "x":
            print("Exit the game.")
            exit(0)
        elif user_input == "n":
            print("Do You want to play another game?")
            break
        elif user_input == "s":
            print(f"The hidden number is {number_to_guess}.")
        else:
            guesses += 1
            user_input = int(user_input)
            if user_input < number_user:
                print("The number is too small!")
            elif user_guess > number_to_guess:
                print("The number is too big!")
            else:
                print(f"Congratulations! You guessed the number in {guesses} attempts.")
                break
def ask_play_again():
    while True:
        play_again = input("\nDo you want to play again? (y/n): ").strip().lower()
        if play_again == "y":
            return True
        elif play_again == "n":
            return False
        else:
            print("Invalid input. Please enter 'y' or 'n'.")

def get_user_input():
    while True:
        user_input = input("Enter your guess (or 'x' to exit, 'n' to start a new game, 's' to show the number): ").strip().lower()
        if user_input in {"x", "n", "s"}:
            return user_input
        elif user_input.isdigit():
            if 1 <= int(user_input) <= 20:
                return user_input
            else:
                print("Invalid number. Please enter a number between 1 and 20.")
        else:
            print("Invalid input. Please try again.")

if __name__ == "__main__":
    main()
