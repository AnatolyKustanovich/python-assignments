import random

def main():
    while True:
        if not game():
            break
    print("Thanks You! Bye!")

def game():
    number = random.randint(1, 20)
    attempts = 0

    print("\nPlease enter a number between 1 and 20!")
    while True:
        you_input = input("Enter your guess (You can use 'x' to exit, 'n' for a new game, 's' to show the number): ").strip().lower()

        if you_input == "x":
            print("Exiting the game. Goodbye!")
            exit(0)
        elif you_input == "n":
            print("You chose to abandon this game.")
            return play_again()
        elif you_input == "s":
            print(f"The hidden number is {number}.")
        elif you_input.isdigit():
            number_user = int(you_input)
            attempts += 1

            if number_user > number:
                print("Your number is too big. Enter another.")
            elif number_user < number:
                print("Your number is too small. Enter another.")
            else:
                print(f"Exactly! The correct number is {number}.")
                print(f"You guessed it in {attempts} attempts.")
                return play_again()
        else:
            print("Invalid input! Please enter a number between 1 and 20, or 'x', 'n', 's'.")

def play_again():
    while True:
        play_again = input("\nDo you want to play a new game? (y/n): ").strip().lower()
        if play_again == "y":
            return True
        elif play_again == "n":
            return False
        else:
            print("Invalid input. Please enter 'y' or 'n'.")

if __name__ == "__main__":
    main()
