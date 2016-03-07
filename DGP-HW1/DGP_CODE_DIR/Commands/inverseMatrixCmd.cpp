#include "stdafx.h"
#include "Utils/MatlabInterface.h"
#include "Utils/Maya_Macros.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/Utilities.h"
#include "Utils/GMM_Macros.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "inverseMatrixCmd.h"

inverseMatrixCmd::inverseMatrixCmd()
{
}

void * inverseMatrixCmd::creator()
{
	return new inverseMatrixCmd;
}

MSyntax inverseMatrixCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	for (int i = 0; i < 9; i++) {//get 9 double variables
		stat = commandSyntax.addArg(MSyntax::kDouble);
		MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	}

	return commandSyntax;
}

MString inverseMatrixCmd::commandName()
{
	return "inverseMatrixCmd";
}

bool inverseMatrixCmd::isUndoable() const
{
	return false;
}

MStatus inverseMatrixCmd::doIt(const MArgList & argList)
{
	MStatus stat = MS::kSuccess;//returned status
	MSyntax commandSyntax = syntax();// chack that the syntax is correct
	MArgParser argData(commandSyntax, argList, &stat);//if the syntax is correct, insert argList to argData
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());// if it is not correct - error message

	GMMDenseColMatrix mToInverse(3,3); //a gmm matrix, dimensions are: 3 x 3

	for (int i = 0; i < argList.length(); i++) {//fill the matrix with the arguments from the command
		int row = i / 3;
		int column = i % 3;
		stat = argData.getCommandArgument(i, mToInverse(row, column));
		MCHECKERROR(stat, "Can't access object list");
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("toInverse", mToInverse);//send the matrix to matlab

	//load the inverse script
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/inverse.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'inverse.m' failed with error code " << res << std::endl;
	}

	std::cerr << "before inverse " << mToInverse << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("toInverse", mToInverse);//get the incersed matrix from matlab
	
	//std::cout <<"after inverse " << mToInverse << std::endl;
	std::cerr << "after inverse " << mToInverse << std::endl;

	cout.flush();

	return stat;
}