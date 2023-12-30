set nocp
set backspace=indent,eol,start
set mouse=a

set autoindent
set shiftwidth=4

set laststatus=2
set statusline=\ %F\ %=Position\ %=%l:%v\

" sudo apt-get install vim-gui-common
" sudo apt-get install vim-runtime
" Highlight the syntax depending on a file format
if has("syntax")
syntax on
endif

" Jump to the last position when reopenming a file
if has("autocmd")
au BufReadPost * if line("'\"") > 0 && line("'\"") <= line("$")
\| exe "normal! g'\"" | endif
endif
