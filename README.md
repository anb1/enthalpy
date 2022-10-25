# enthalpy
Calculate enthalpy of a target reaction and set of additional reactions.


Input

	"Please enter target reaction:"
	User: 	P4O10 + 6H2O -> 4H3PO4 
	"Please enter additional reactions, followed by enthalpy in kilojoules, 
		Use space and semicolon to separate reaction and enthalpy, 
		Ctrl-Z in Windows/Ctrl-D in Unix to end input:"
	User: 	P4 + 5O2 -> P4O10 ; -2984
		2H2 + O2 -> 2H2O ; -571.6
		P4 + 6H2 + 8O2 -> 4H3PO4 ; -2568.8
		(ctrl-z/ctrl-d)


Process

	(1) Target reaction 
    Overview: "P4O10 + 6H2O -> 4H3PO4" should be interpreted as the linear expression x_1+6x_2-4x_3.
			4H3PO4 has a negative coefficient because it's to the right of the arrow.
		Separate the coefficients from the chemical formulae (using split() and stoi()).
			Put chemical formulae in vector<string> formula.
				formula == {"P4O10", "H2O", "H3PO4"}
			Put coefficients in vector<int> targetRxnCoefficients.
				targetRxnCoefficients == {1, 6, -4}
    Target reaction gets its own vector for reasons that will become clear later.
	(2) Additional reactions 
    Overview: "P4 + 5O2 -> P4O10 ; -2984" should be interpreted as the linear equation x_4+ 5x_5– x_1= -2984.
      (Just like the target reaction, coefficients to the right of the arrow should be negative.)
		Separate the coefficients from the chemical formulae (using split() and stoi()).
			Append chemical formulae in vector<string> formula. 
				formula == {"P4O10", "H2O", "H3PO4", "P4", "O2"}
			Put coefficients in vector<vector<int>> reactions. 
				reactions is a 2-dimensional vector. 
          The row number corresponds to the reaction number (e.g., rxn1, rxn2, rxn3, etc.), minus one, to index.
          The column number corresponds to the position of the compound in vector<string> formula.
        "P4 + 5O2 -> P4O10 ; -2984" (rxn1) is interpreted to mean 
        reactions[0] == {-1, 0, 0, 1, 5}.
          The zeroes indicate H2O and H3PO4, which are not present but maintain their positions.
			Store enthalpy in vector<double> b.
        The position in b also corresponds to the reaction number (minus one, to index).
				"P4 + 5O2 -> P4O10 ; -2984" (rxn1) is interpreted to mean b[0] == -2984.
		Arrange all the coefficients together in a matrix A of dimensions n x m. 
			Matrix construction
				n == targetRxnCoefficients.size() == 3, in this example.
				m == 1 + the final formula.size().
          formula.size() is the number of x vars or unique chemical compounds from all the reactions.
          Plus 1, because this is an augmented matrix and the last column should be the enthalpies (dH) in vector b.
			The first row of matrix A is {-1, 0, 0, 1, 5, 0, -2894}.
        This is the same as reactions[0] earlier except for the additional 0 in column 6, and dH at the end.
          (Zeroes should fill all columns between the final column of reactions[i] and dH.)
		When matrix A is finally constructed, the matrix should be row reduced.
			In this example, the first 3 rows and columns of the matrix should look like:
				{1, 0, 0, ... }
				{0, 1, 0, ... }
				{0, 0, 1, ... }
				This is row-reduced echelon form. 
    Each row in the row-reduced matrix should be multiplied by the corresponding entry in targetRxnCoefficients.
			targetRxnCoefficients == {1, 6, -4}; therefore, matrix A will look like this: 
				{1, 0, 0, 1*(...) }
				{0, 6, 0, 6*(...) }
				{0, 0, -4, -4*(...) }
		All the rows of matrix A should be summed.
			E.g. {1, 6, -4, 1*(...) + 6*(...) + -4*(...) }
			
			
Output
	Print to standard output everything in matrix A after column n.
		E.g. print the portion that is 1*(...) + 6*(...) + -4*(...) from the example above.
	If the chemical reactions provided to the program are consistent, then all the variables should cancel out and only scalar values of enthalpy (dH) should be left.
		E.g. 2130 kilojoules
