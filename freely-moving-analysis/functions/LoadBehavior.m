function [ BehaviorFiles ] = LoadBehavior(varargin)
    
%LoadBehavior(BehaviorFilePath, AnimalList, DateRange, SessionNumber, Variables)
%This function uses the BehaviorFilePath, AnimalList, DateRange, and SessionNumber inputs to filter one or more .Mtrx  behavior data files. It will 
%load the relevant files and check if their trial parameters match the other
%input arguments. If they do, then it will return the loaded file.

%Inputs:
%
%   BehaviorFilePath: path to directory with behavior files
%
%   AnimalList: String of animal IDs separated by commas and no spaces
%   (e.g., 'SZ001,SZ002,SZ003') that are used to filter behavior files
%
%   DateRange: array of dates (e.g., [YYYYMMDD, ..., YYYYMMDD]) that are used to filter behavior files
%
%   SessionNumber: array of session numbers (e.g., [1]) that are used to filter behavior files
%
%   Variables: comma delimited string of variables (e.g.,
%   'puffData,config,experiment') that the function will attempt to load
%   from each filtered behavior files. If a Variables input is not included,
%   all variables will be loaded.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:numel(varargin)
        if strcmp(varargin(i),{'BehaviorFilePath'})
            BehaviorFilePath=cell2mat(varargin(i+1));
        end
        if strcmp(varargin(i),{'AnimalList'})
            AnimalList=strsplit(cell2mat(varargin(i+1)),',');
        end
        if strcmp(varargin(i),{'DateRange'})
            DateRange=cell2mat(varargin(i+1));
        end
        if strcmp(varargin(i),{'SessionNumber'})
            SessionNumber=cell2mat(varargin(i+1));
        end
        if strcmp(varargin(i),{'Variables'})
            Variables=strsplit(cell2mat(varargin(i+1)),',');
        end
    end
    
    % if Variables input was not included, set loadAllVariablesFlag to true
    if exist('Variables','var')
        loadAllVariablesFlag = false;
    else
        loadAllVariablesFlag = true;
    end
    
    BehaviorFileNames = subdir(fullfile(BehaviorFilePath,'*Mtrx.mat'));
    BehaviorFileNamesTrim=[];
    StructNum=1;
    
    if isempty(BehaviorFileNames)
        error('No Mtrx files in this directory');
    end
    
    
 %Identify session files for the correct Animal/Date
    for i=1:length(BehaviorFileNames)
        [~,ShortBehaviorFileName,~] = fileparts(BehaviorFileNames(i).name);
        for j=1:length(AnimalList)
            if contains(ShortBehaviorFileName,AnimalList(j))
                if ~isempty(find(DateRange == str2double(ShortBehaviorFileName(1:8)),1))
                    if ~isempty(SessionNumber)
                        for k =1:length(SessionNumber)
                            if ~isempty(find(SessionNumber(k) == str2double(ShortBehaviorFileName(10:11)),1))
                                BehaviorFileNamesTrim=[BehaviorFileNamesTrim, BehaviorFileNames(i,:)];
                            end
                        end
                    else
                        BehaviorFileNamesTrim=[BehaviorFileNamesTrim, BehaviorFileNames(i,:)];
                    end
                end
            end
        end
    end
    
    
    %Load files and make sure they match the requested parameter options
    
    if isempty(BehaviorFileNamesTrim)
        error('No Mtrx files matching those criteria');
    end
    
    for i=1:length(BehaviorFileNamesTrim)
        
        % Check if all variables should be loaded
        if loadAllVariablesFlag
            Variables = who('-file', BehaviorFileNamesTrim(i).name);
        end

        for j = 1:length(Variables)
            load(BehaviorFileNamesTrim(i).name, Variables{j});
            BehaviorFiles(i).name = BehaviorFileNamesTrim(i).name;
            BehaviorFiles(i).(Variables{j}) = eval(Variables{j});
        end
        fprintf('\n %s \n', BehaviorFileNamesTrim(i).name)
     
    if ~exist('BehaviorFiles','var')
        fprintf('\nThere are no behavior files matching your parameters\n');
        BehaviorFiles=[];
    end

    end

