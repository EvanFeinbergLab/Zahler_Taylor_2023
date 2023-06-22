function [ BehaviorFiles ] = LoadBehavior(varargin)
    
%LoadBehavior(BehaviorFilePath, AnimalList, DateRange, SessionLengthRange, ITIRange, MinAngleRange, Quiescence, BiasCorrection, VisualFeedback, AuditoryFeedback)
%This function will take as input the name of a behavior file and a vector of dates. It will 
%load the relevant files and check if their trial parameters match the other
%input arguments. If they do, then it will return the loaded file.

%Input:
    %BehaviorFilePath=name of location of behavior directory.
    
    %AnimalList= String of animal IDs separated by commas and no spaces ex. 'SZ001,SZ002,SZ003'
    
    %DateRange= array of dates
    
    %Comments
    
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
    
    BehaviorFileNames = subdir(fullfile(BehaviorFilePath,'*Mtrx.mat'));
    BehaviorFileNamesTrim=[];
    StructNum=1;
    
    
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
    
    for i=1:length(BehaviorFileNamesTrim)

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

