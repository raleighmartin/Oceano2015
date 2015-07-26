%% function to determine the type (i.e. string, double, etc.) of variables in 'tableVariableList'
% based on information about all variables in 'InstrumentVariables'
% function outputs a list of types corresponding to the input elements in 'tableVariableList'
% Dependencies: NONE
% Used by: LoggerImport

function tableDataTypes = GetDataTypes(tableVariableList,InstrumentVariables)

%how many variables?
N_variables = length(tableVariableList);

%inialize list of types
tableDataTypes = cell(N_variables,1);

%go through each variable
for i = 1:N_variables
    Variable = tableVariableList{i}; %get variable name
    tableIndex = find(strcmp(InstrumentVariables.VarNameSpecific,Variable));
    %tableIndex = CellStrFind(InstrumentVariables.VarNameSpecific,Variable); %get index of Variable in metadata file
    tableDataTypes{i} = char(InstrumentVariables.DataType(tableIndex)); %get type for this index
end