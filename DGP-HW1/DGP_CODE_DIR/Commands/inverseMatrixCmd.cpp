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

	//create a syntax that takes 9 doubles
	for (int i = 0; i < 9; i++) {
		stat = commandSyntax.addArg(MSyntax::kDouble);
		MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	}

	return commandSyntax;
}

MString inverseMatrixCmd::commandName()
{
	return "inverseMatrixCmd";
}

bool inverseMatrixCmd::isUndoable(const MArgList & argList) const
{
	MStatus stat = MS::kSuccess;//returned status
	MSyntax commandSyntax = syntax();// chack that the syntax is correct
	MArgParser argData(commandSyntax, argList, &stat);//if the syntax is correct, insert argList to argData
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());// if it is not correct - error message

	GMMDenseColMatrix matToInverse(3, 3); //a gmm matrix, dimensions are: 3 x 3
	GMMDenseColMatrix isInvertable(1, 1); //a gmm matrix, dimensions are: 1 x 1

										  //fill the matrix with the arguments from the command
	for (int i = 0; i < argList.length(); i++) {
		int row = i / 3;
		int column = i % 3;
		stat = argData.getCommandArgument(i, matToInverse(row, column));
		MCHECKERROR(stat, "Can't access object list");
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("matToInverse", matToInverse);//send the matrix to matlab

	//load the issingular script
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/DGP/DGP-Ben/matlab scripts/isSingular.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'inverse.m' failed with error code " << res << std::endl;
	}

	std::cerr << "before singularity check " << matToInverse << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("isInvertable", isInvertable);//get the incersed matrix from matlab

	std::cerr << "after singularity check " << isInvertable << std::endl;
	cout.flush();

	if (isInvertable(0,0)==0)
	{
		return true;
	}
	else 
	{
		return false;
	}
}

MStatus inverseMatrixCmd::doIt(const MArgList & argList)
{
	MStatus stat = MS::kSuccess;//returned status
	MSyntax commandSyntax = syntax();// chack that the syntax is correct
	MArgParser argData(commandSyntax, argList, &stat);//if the syntax is correct, insert argList to argData
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());// if it is not correct - error message

	if (isUndoable(argList)) {
		stat = MS::kFailure;
		MCHECKERROR(stat, "Matrix is sungular");
	}

	GMMDenseColMatrix matToInverse(3,3); //a gmm matrix, dimensions are: 3 x 3
	GMMDenseColMatrix inveresedMat(3, 3); //a gmm matrix, dimensions are: 3 x 3

	//fill the matrix with the arguments from the command
	for (int i = 0; i < argList.length(); i++) {
		int row = i / 3;
		int column = i % 3;
		stat = argData.getCommandArgument(i, matToInverse(row, column));
		MCHECKERROR(stat, "Can't access object list");
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("matToInverse", matToInverse);//send the matrix to matlab

	//load the inverse script
	//int res = MatlabInterface::GetEngine().LoadAndRunScript("%DGP_CODE_DIR%/matlab scripts/inverse.m");
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/DGP/DGP-Ben/matlab scripts/inverse.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'inverse.m' failed with error code " << res << std::endl;
	}

	std::cerr << "before inverse " << matToInverse << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("inversedMat", inveresedMat);//get the incersed matrix from matlab
	
	//std::cout <<"after inverse " << mToInverse << std::endl;
	std::cerr << "after inverse " << inveresedMat << std::endl;

	cout.flush();

	return stat;
}