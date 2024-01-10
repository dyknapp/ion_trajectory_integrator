% Check if there are large files that can't be saved to git.

filelist = dir(fullfile('.', '**\*.*'));  %get list of files and directories in any subfolder
filelist = filelist(~[filelist.isdir]); % remove directories from list

git_storage_size_limit = 1 * (2^20);

found_ignores = false;
ignore_files = [" "];
ignores = 1;
for idx = 1:length(filelist)
    file = filelist(idx);
    path = extractAfter(strcat(file.folder, "\", file.name), pwd);

    if ~contains(path, '.git')
        size = dir("." + path).bytes;
        fprintf('%d bytes:\n%s\n', size, path);
        if size >= git_storage_size_limit
            fprintf("Ignoring this file for git because it is too large.\n")
            ignore_files(ignores) = strrep(extractAfter(path, "\"), "\", "/");
            ignores = ignores + 1;
            found_ignores = true;
        end
    end
end

% Update gitignore file
fileID = fopen('.gitignore', 'r');
gitignore_boilerplate_text = [extractBefore(fscanf(fileID, '%c'), "# END OF BOILERPLATE"), newline, '# END OF BOILERPLATE', newline, '# LARGE FILES:'];
fclose(fileID);
for ignore_addition = ignore_files
    gitignore_boilerplate_text = string(gitignore_boilerplate_text) + string(newline + ignore_addition);
end
fileID = fopen('.gitignore', 'w');
fprintf(fileID, '%s', gitignore_boilerplate_text);
fclose(fileID);

fprintf("\nIgnored %d files.\n", ignores - 1);