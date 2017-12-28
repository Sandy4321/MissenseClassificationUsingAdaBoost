function err = calculate_error(truey, predy)
	size(truey)
	size(predy)
	err = sum(truey==predy)/size(truey,1);
end
